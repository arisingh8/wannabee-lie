/*
 * This Kotlin source file was generated by the Gradle 'init' task.
 */
package wannabee.lie

import io.kotest.core.spec.style.StringSpec
import io.kotest.property.Arb
import io.kotest.property.arbitrary.bind
import io.kotest.property.arbitrary.double
import io.kotest.property.arbitrary.map
import io.kotest.property.forAll
import org.ejml.data.DMatrix2
import org.ejml.data.DMatrix
import org.ejml.data.DMatrix3
import org.ejml.data.DMatrix3x3
import org.ejml.dense.fixed.CommonOps_DDF3
import kotlin.math.PI

fun DMatrix.equalsDelta(other: DMatrix): Boolean {
    if (this.numCols != other.numCols || this.numRows != other.numRows) {
        return false
    }
    for (row in 0..<this.numRows) {
        for (col in 0..<this.numCols) {
            if (!this.get(row, col).equalsDelta(other.get(row, col))) {
                return false
            }
        }
    }
    return true
}

class GeometryTest: StringSpec({
    val rotationArb: Arb<LieRotation2d> = Arb.double(-PI..PI)
        .map { LieRotation2d.exp(it) }
    val vectorArb: Arb<DMatrix2> = Arb.bind(
        Arb.double(-1000.0..1000.0),
        Arb.double(-1000.0..1000.0)
    ) { x, y -> DMatrix2(x, y) }
    val twistArb: Arb<LieTwist2d> = Arb.bind(
        vectorArb,
        Arb.double(-1000.0..1000.0)
    ) { vector, angle -> LieTwist2d(vector, angle) }
    val poseArb: Arb<LiePose2d> = Arb.bind(
        vectorArb,
        rotationArb
    ) { vector, heading -> LiePose2d(vector, heading) }

    "rotation composed with the identity should equal itself" {
        forAll(rotationArb) { rot ->
            (rot * LieRotation2d()).equalsDelta(rot)
        }
    }

    "rotation composed with its inverse should equal the identity rotation" {
        forAll(rotationArb) { rotation ->
            val inverse = rotation.inverse()
            (rotation * inverse).equalsDelta(LieRotation2d())
        }
    }

    "rotation composition should be associative" {
        forAll(rotationArb, rotationArb, rotationArb) { rot1, rot2, rot3 ->
            ((rot1 * rot2) * rot3).equalsDelta(rot1 * (rot2 * rot3))
        }
    }

    "rotation identity acting on points should yield the same points" {
        forAll(vectorArb) { vector ->
            (LieRotation2d() * vector).equalsDelta(vector)
        }
    }

    "rotation group actions should be compatible" {
        forAll(rotationArb, rotationArb, vectorArb) { rot1, rot2, vector ->
            ((rot1 * rot2) * vector).equalsDelta(rot1 * (rot2 * vector))
        }
    }

    "rotation exp and log should invert each other" {
        forAll(Arb.double(-PI..PI)) { angle ->
            LieRotation2d.exp(angle).log().equalsDelta(angle)
        }
    }

    "pose composed with the identity should equal itself" {
        forAll(poseArb) { pose ->
            (pose * LiePose2d()).equalsDelta(pose)
        }
    }

    "pose composed with its inverse should equal the identity pose" {
        forAll(poseArb) { pose ->
            (pose * pose.inverse()).equalsDelta(LiePose2d())
        }
    }

    "pose composition should be associative" {
        forAll(poseArb, poseArb, poseArb) { pose1, pose2, pose3 ->
            ((pose1 * pose2) * pose3).equalsDelta(pose1 * (pose2 * pose3))
        }
    }

    "pose identity acting on points should yield the same points" {
        forAll(vectorArb) { vector ->
            (LiePose2d() * vector).equalsDelta(vector)
        }
    }

    "pose group actions should be compatible" {
        forAll(poseArb, poseArb, vectorArb) { pose1, pose2, vector ->
            ((pose1 * pose2) * vector).equalsDelta(pose1 * (pose2 * vector))
        }
    }

    "exp satisfies combinations" {
        forAll(twistArb, Arb.double(-1000.0, 1000.0), Arb.double(-1000.0, 1000.0)) { twist, t, s ->
            val tauVec = DMatrix3(twist.translation.a1, twist.translation.a2, twist.angle)
            val tau1 = DMatrix3()
            val tau2 = DMatrix3()
            val tau3 = DMatrix3()
            CommonOps_DDF3.scale(t + s, tauVec, tau1)
            CommonOps_DDF3.scale(t, tauVec, tau2)
            CommonOps_DDF3.scale(s, tauVec, tau3)
            LiePose2d.exp(LieTwist2d(tau1.a1, tau1.a2, tau1.a3)).equalsDelta(LiePose2d.exp(LieTwist2d(tau2.a1, tau2.a2, tau2.a3)) * LiePose2d.exp(LieTwist2d(tau3.a1, tau3.a2, tau3.a3)))
        }
    }

    "exp inverse property" {
        forAll(twistArb) { twist ->
            val tauVec = DMatrix3(twist.translation.a1, twist.translation.a2, twist.angle)
            val tau = DMatrix3()
            CommonOps_DDF3.scale(-1.0, tauVec, tau)
            LiePose2d.exp(LieTwist2d(tau.a1, tau.a2, tau.a3)).equalsDelta(LiePose2d.exp(LieTwist2d(tauVec.a1, tauVec.a2, tauVec.a3)).inverse())
        }
    }

    "adjoint plus thingy should be consistent" {
        forAll(poseArb, twistArb) { pose, twist ->
            val tauVec = DMatrix3(twist.translation.a1, twist.translation.a2, twist.angle)
            val adjProd = DMatrix3()
            CommonOps_DDF3.mult(pose.adjoint(), tauVec, adjProd)
            val newTwist = LieTwist2d(adjProd.a1, adjProd.a2, adjProd.a3)
            (pose + twist).equalsDelta(LiePose2d.exp(newTwist) * pose)
        }
    }

    "adjoint inverse is consistent" {
        forAll(poseArb) { pose ->
            val inverseAdj = DMatrix3x3()
            val invert = CommonOps_DDF3.invert(pose.adjoint(), inverseAdj)
            if (!invert) return@forAll true
            pose.inverse().adjoint().equalsDelta(inverseAdj)
        }
    }

    "adjoint composition is consistent" {
        forAll(poseArb, poseArb) { pose1, pose2 ->
            val adjProd = DMatrix3x3()
            CommonOps_DDF3.mult(pose1.adjoint(), pose2.adjoint(), adjProd)
            (pose1 * pose2).adjoint().equalsDelta(adjProd)
        }
    }

    "inverse Jacobian matches negative adjoint" {
        forAll(poseArb) { pose ->
            val (_, jacobian) = pose.inverseJacobians()
            val negativeAdjoint = pose.adjoint()
            CommonOps_DDF3.changeSign(negativeAdjoint)
            jacobian.equalsDelta(negativeAdjoint)
        }
    }

    "composition Jacobian matches adjoint of other" {
        forAll(poseArb, poseArb) { pose1, pose2 ->
            val (_, jacobian) = pose1.timesJacobians(pose2)
            val invAdjoint = DMatrix3x3()
            val invert = CommonOps_DDF3.invert(pose2.adjoint(), invAdjoint)
            if (!invert) return@forAll true
            jacobian.equalsDelta(invAdjoint)
        }
    }

    "right and left Jacobians match each other" {
        forAll(twistArb) { twist ->
            LieTwist2d(-twist.translation.a1, -twist.translation.a2, -twist.angle).rightJ().equalsDelta(twist.leftJ())
        }
    }

    "plus Jacobians match adjoint and right Jacobian" {
        forAll(poseArb, twistArb) { pose, twist ->
            val (_, selfJ, otherJ) = pose.plusJacobians(twist)
            val invAdjoint = DMatrix3x3()
            CommonOps_DDF3.invert(LiePose2d.exp(twist).adjoint(), invAdjoint)
            selfJ.equalsDelta(invAdjoint) && otherJ.equalsDelta(twist.rightJ())
        }
    }
})