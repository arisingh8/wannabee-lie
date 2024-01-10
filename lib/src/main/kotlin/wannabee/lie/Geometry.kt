package wannabee.lie

import org.ejml.data.DMatrix2
import org.ejml.data.DMatrix2x2
import org.ejml.data.DMatrix3x3
import org.ejml.data.DMatrixRMaj
import org.ejml.dense.fixed.CommonOps_DDF2
import org.ejml.dense.fixed.CommonOps_DDF3
import org.ejml.dense.fixed.NormOps_DDF2
import kotlin.math.*

@JvmOverloads
fun skew(scalar: Double = 1.0) = DMatrix2x2(
    0.0, -1*scalar,
    scalar, 0.0
)

class LieVector2d(val vector: DMatrix2 = DMatrix2()) {
    val x get() = vector.a1
    val y get() = vector.a2

    constructor(x: Double, y: Double) : this(DMatrix2(x, y))

    fun norm() = NormOps_DDF2.fastNormF(vector)
    operator fun component1() = x
    operator fun component2() = y
    override fun toString() = "LieVector2d(x: $x y: $y)"
}

class LieRotation2d(val rotMatrix: DMatrix2x2) {
    companion object {
        @JvmStatic
        fun exp(heading: Double): LieRotation2d {
            return LieRotation2d(DMatrix2x2(
                cos(heading), -sin(heading),
                sin(heading), cos(heading)
            ))
        }
    }

    fun log() = atan2(rotMatrix.a21, rotMatrix.a11)
    fun inverse(): LieRotation2d {
        val tmp = DMatrix2x2()
        CommonOps_DDF2.transpose(rotMatrix, tmp)
        return LieRotation2d(tmp)
    }
    fun compose(other: LieRotation2d): LieRotation2d {
        val tmp = DMatrix2x2()
        CommonOps_DDF2.mult(rotMatrix, other.rotMatrix, tmp)
        return LieRotation2d(tmp)
    }
    fun compose(other: LieVector2d): LieVector2d {
        val tmp = DMatrix2()
        CommonOps_DDF2.mult(rotMatrix, other.vector, tmp)
        return LieVector2d(tmp)
    }
}

data class LieTwist2d(val translation: LieVector2d, val angle: Double) {
    constructor(x: Double, y: Double, heading: Double) : this(LieVector2d(x, y), heading)

    fun rightJ(): DMatrix3x3 {
        return if (angle == 0.0) {
            DMatrix3x3(
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 1.0
            )
        } else {
            DMatrix3x3(
                sin(angle), (1 - cos(angle))/angle, (angle*translation.x - translation.y + translation.y*cos(angle) - translation.x*sin(angle))/angle.pow(2),
                (cos(angle) - 1)/angle, sin(angle), (translation.x + angle*translation.y - translation.x*cos(angle) - translation.y*sin(angle))/angle.pow(2),
                0.0, 0.0, 1.0
            )
        }
    }

    fun leftJ(): DMatrix3x3 {
        return if (angle == 0.0) {
            DMatrix3x3(
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 1.0
            )
        } else {
            DMatrix3x3(
                sin(angle), (cos(angle) - 1)/angle, (angle*translation.x + translation.y - translation.y*cos(angle) - translation.x*sin(angle))/angle.pow(2),
                (1 - cos(angle))/angle, sin(angle), (-translation.x + angle*translation.y + translation.x*cos(angle) - translation.y*sin(angle))/angle.pow(2),
                0.0, 0.0, 1.0
            )
        }
    }
}

data class LiePose2d(val position: LieVector2d, val rotation: LieRotation2d) {
    constructor(position: LieVector2d, heading: Double) : this(position, LieRotation2d.exp(heading))
    constructor(x: Double, y: Double, heading: Double) : this(LieVector2d(x, y), heading)

    companion object {
        @JvmStatic
        fun exp(twist: LieTwist2d): LiePose2d {
            val v = DMatrix2x2(
                sin(twist.angle), -(1 - cos(twist.angle)),
                1 - cos(twist.angle), sin(twist.angle)
            )
            if (twist.angle != 0.0) {
                CommonOps_DDF2.scale(1 / twist.angle, v)
            }
            val position = DMatrix2()
            CommonOps_DDF2.mult(v, twist.translation.vector, position)
            return LiePose2d(LieVector2d(position), LieRotation2d.exp(twist.angle))
        }
    }

    fun log(): LieTwist2d {
        val angle = rotation.log()
        val vinv = DMatrix2x2()
        if (angle == 0.0) {
            CommonOps_DDF2.setIdentity(vinv)
        } else {
            val a = sin(angle) / angle
            val b = (1 - cos(angle)) / angle
            vinv.setTo(
                a, b,
                -b, a
            )
            CommonOps_DDF2.scale(1/(a.pow(2) + b.pow(2)), vinv)
        }
        val translation = DMatrix2()
        CommonOps_DDF2.mult(vinv, position.vector, translation)
        return LieTwist2d(LieVector2d(translation), angle)
    }

    fun adjoint(): DMatrix3x3 {
        val tmp = DMatrix2()
        CommonOps_DDF2.mult(skew(), position.vector, tmp)
        CommonOps_DDF2.changeSign(tmp)
        return DMatrix3x3(
            rotation.rotMatrix.a11, rotation.rotMatrix.a12, tmp.a1,
            rotation.rotMatrix.a21, rotation.rotMatrix.a22, tmp.a2,
            0.0, 0.0, 1.0
        )
    }

    fun plus(other: LieTwist2d): PlusResult {
        val tau = LiePose2d.exp(other)
        val result = compose(tau)
        return PlusResult(
            result.pose,
            tau.inverse().pose.adjoint(),
            other.rightJ()
        )
    }

    fun inverse(): InverseResult {
        val inverseRotation = rotation.inverse()
        val inverseTranslation = inverseRotation.compose(position)
        CommonOps_DDF2.changeSign(inverseTranslation.vector)

        val jSelf = adjoint()
        CommonOps_DDF3.changeSign(jSelf)
        return InverseResult(
            LiePose2d(inverseTranslation, inverseRotation),
            jSelf
        )
    }

    fun compose(other: LieVector2d): PointsComposeResult {
        val result = DMatrix2()
        CommonOps_DDF2.mult(rotation.rotMatrix, other.vector, result)
        CommonOps_DDF2.addEquals(result, position.vector)

        val tmp = DMatrix2()
        val tmp2 = DMatrix2x2()
        CommonOps_DDF2.mult(rotation.rotMatrix, skew(), tmp2)
        CommonOps_DDF2.mult(tmp2, other.vector, tmp)
        val jSelf = DMatrixRMaj(arrayOf(
            doubleArrayOf(rotation.rotMatrix.a11, rotation.rotMatrix.a12, tmp.a1),
            doubleArrayOf(rotation.rotMatrix.a21, rotation.rotMatrix.a22, tmp.a2)
        ))
        return PointsComposeResult(
            LieVector2d(result),
            jSelf,
            rotation.rotMatrix.copy()
        )
    }
    fun compose(other: LiePose2d): PoseComposeResult {
        val newTranslation = compose(other.position).point
        val newRotation = rotation.compose(other.rotation)
        return PoseComposeResult(
            LiePose2d(newTranslation, newRotation),
            other.inverse().pose.adjoint()
        )
    }
}

data class InverseResult(val pose: LiePose2d, val jSelf: DMatrix3x3)
data class PlusResult(val pose: LiePose2d, val jSelf: DMatrix3x3, val jTau: DMatrix3x3)
data class PointsComposeResult(val point: LieVector2d, val jSelf: DMatrixRMaj, val jTau: DMatrix2x2)
data class PoseComposeResult(val pose: LiePose2d, val jSelf: DMatrix3x3)