package wannabee.lie

import org.ejml.data.DMatrix2
import org.ejml.data.DMatrix2x2
import org.ejml.data.DMatrix3
import org.ejml.data.DMatrix3x3
import org.ejml.data.DMatrixRMaj
import org.ejml.dense.fixed.CommonOps_DDF2
import org.ejml.dense.fixed.CommonOps_DDF3
import org.ejml.dense.fixed.NormOps_DDF2
import kotlin.math.*

private fun skew(scalar: Double = 1.0) = DMatrix2x2(
    0.0, -1*scalar,
    scalar, 0.0
)

fun Double.equalsDelta(other: Double) = abs(this - other) < 0.00000001
fun DMatrix2.equalsDelta(other: DMatrix2) = this.a1.equalsDelta(other.a1) && this.a2.equalsDelta(other.a2)

class LieRotation2d(
    val rotMatrix: DMatrix2x2 = DMatrix2x2(
        1.0, 0.0,
        0.0, 1.0
    )
) {
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
    operator fun times(other: LieRotation2d): LieRotation2d {
        val tmp = DMatrix2x2()
        CommonOps_DDF2.mult(rotMatrix, other.rotMatrix, tmp)
        return LieRotation2d(tmp)
    }
    operator fun times(other: DMatrix2): DMatrix2 {
        val result = DMatrix2()
        CommonOps_DDF2.mult(rotMatrix, other, result)
        return result
    }

    override fun toString() = "LieRotation2d(angle=${log()})"

    fun equalsDelta(other: LieRotation2d) = log().equalsDelta(other.log())
}

data class LieTwist2d(
    val translation: DMatrix2 = DMatrix2(0.0, 0.0),
    val angle: Double = 0.0
) {
    constructor(x: Double, y: Double, heading: Double) : this(DMatrix2(x, y), heading)

    fun rightJ(): DMatrix3x3 {
        return if (angle == 0.0) {
            DMatrix3x3(
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 1.0
            )
        } else {
            DMatrix3x3(
                sin(angle)/angle, (1 - cos(angle))/angle, (angle*translation.a1 - translation.a2 + translation.a2*cos(angle) - translation.a1*sin(angle))/angle.pow(2),
                (cos(angle) - 1)/angle, sin(angle)/angle, (translation.a1 + angle*translation.a2 - translation.a1*cos(angle) - translation.a2*sin(angle))/angle.pow(2),
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
                sin(angle)/angle, (cos(angle) - 1)/angle, (angle*translation.a1 + translation.a2 - translation.a2*cos(angle) - translation.a1*sin(angle))/angle.pow(2),
                (1 - cos(angle))/angle, sin(angle)/angle, (-translation.a1 + angle*translation.a2 + translation.a1*cos(angle) - translation.a2*sin(angle))/angle.pow(2),
                0.0, 0.0, 1.0
            )
        }
    }

    override fun toString() = "LieTwist2d(x: ${translation.a1}, y: ${translation.a2}, angle: ${angle})"

    fun equalsDelta(other: LieTwist2d) = translation.equalsDelta(other.translation) && angle.equalsDelta(other.angle)
}

data class LiePose2d(
    val position: DMatrix2 = DMatrix2(0.0, 0.0),
    val rotation: LieRotation2d = LieRotation2d()
) {
    constructor(position: DMatrix2, heading: Double) : this(position, LieRotation2d.exp(heading))
    constructor(x: Double, y: Double, heading: Double) : this(DMatrix2(x, y), heading)

    companion object {
        @JvmStatic
        fun exp(twist: LieTwist2d): LiePose2d {
            if (twist.angle == 0.0) {
                return LiePose2d(twist.translation.copy(), LieRotation2d.exp(twist.angle))
            }
            val v = DMatrix2x2(
                sin(twist.angle) / twist.angle, -(1 - cos(twist.angle)) / twist.angle,
                (1 - cos(twist.angle)) / twist.angle, sin(twist.angle) / twist.angle
            )
            val position = DMatrix2()
            CommonOps_DDF2.mult(v, twist.translation, position)
            return LiePose2d(position, LieRotation2d.exp(twist.angle))
        }
    }

    fun log(): LieTwist2d {
        val angle = rotation.log()
        if (angle == 0.0) {
            return LieTwist2d(position, angle)
        }
        val vinv = DMatrix2x2()
        val a = sin(angle) / angle
        val b = (1 - cos(angle)) / angle
        vinv.setTo(
            a, b,
            -b, a
        )
        CommonOps_DDF2.scale(1/(a.pow(2) + b.pow(2)), vinv)

        val translation = DMatrix2()
        CommonOps_DDF2.mult(vinv, position, translation)
        return LieTwist2d(translation, angle)
    }

    fun adjoint(): DMatrix3x3 {
        val tmp = DMatrix2()
        CommonOps_DDF2.mult(skew(), position, tmp)
        CommonOps_DDF2.changeSign(tmp)
        return DMatrix3x3(
            rotation.rotMatrix.a11, rotation.rotMatrix.a12, tmp.a1,
            rotation.rotMatrix.a21, rotation.rotMatrix.a22, tmp.a2,
            0.0, 0.0, 1.0
        )
    }

    operator fun plus(other: LieTwist2d): LiePose2d {
        return times(exp(other))
    }
    fun plusJacobians(other: LieTwist2d): PlusResult {
        val tau = exp(other)
        return PlusResult(
            times(tau),
            tau.inverse().adjoint(),
            other.rightJ()
        )
    }

    fun inverse(): LiePose2d {
        val inverseRotation = rotation.inverse()
        val inverseTranslation = inverseRotation * position
        CommonOps_DDF2.changeSign(inverseTranslation)
        return LiePose2d(inverseTranslation, inverseRotation)
    }
    fun inverseJacobians(): InverseResult {
        val jSelf = adjoint()
        CommonOps_DDF3.changeSign(jSelf)
        return InverseResult(
            inverse(),
            jSelf
        )
    }

    operator fun times(other: DMatrix2): DMatrix2 {
        val result = DMatrix2()
        CommonOps_DDF2.mult(rotation.rotMatrix, other, result)
        CommonOps_DDF2.addEquals(result, position)
        return result
    }
    fun timesJacobians(other: DMatrix2): PointsComposeResult {
        val tmp = DMatrix2()
        val tmp2 = DMatrix2x2()
        CommonOps_DDF2.mult(rotation.rotMatrix, skew(), tmp2)
        CommonOps_DDF2.mult(tmp2, other, tmp)
        val jSelf = DMatrixRMaj(arrayOf(
            doubleArrayOf(rotation.rotMatrix.a11, rotation.rotMatrix.a12, tmp.a1),
            doubleArrayOf(rotation.rotMatrix.a21, rotation.rotMatrix.a22, tmp.a2)
        ))
        return PointsComposeResult(
            times(other),
            jSelf,
            rotation.rotMatrix.copy()
        )
    }

    operator fun times(other: LiePose2d): LiePose2d {
        val newTranslation = times(other.position)
        val newRotation = rotation * other.rotation
        return LiePose2d(newTranslation, newRotation)
    }
    fun timesJacobians(other: LiePose2d): PoseComposeResult {
        return PoseComposeResult(
            times(other),
            other.inverse().adjoint()
        )
    }

    override fun toString() = "LiePose2d(x: ${position.a1}, y: ${position.a2}, angle: ${rotation})"

    fun equalsDelta(other: LiePose2d) = position.equalsDelta(other.position) && rotation.equalsDelta(other.rotation)
}

data class InverseResult(val pose: LiePose2d, val jSelf: DMatrix3x3)
data class PlusResult(val pose: LiePose2d, val jSelf: DMatrix3x3, val jTau: DMatrix3x3)
data class PointsComposeResult(val point: DMatrix2, val jSelf: DMatrixRMaj, val jTau: DMatrix2x2)
data class PoseComposeResult(val pose: LiePose2d, val jSelf: DMatrix3x3)