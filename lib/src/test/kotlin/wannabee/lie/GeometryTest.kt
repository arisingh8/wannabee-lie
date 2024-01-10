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
import kotlin.math.PI

class GeometryTest: StringSpec({
    val rotationArb: Arb<LieRotation2d> = Arb.double(-PI..PI)
        .map { LieRotation2d.exp(it) }
    val vectorArb: Arb<DMatrix2> = Arb.bind(
        Arb.double(-1000.0..1000.0),
        Arb.double(-1000.0..1000.0)
    ) { x, y -> DMatrix2(x, y) }
    val twistArb: Arb<LieTwist2d> = Arb.bind(
        vectorArb,
        Arb.double(-PI..PI)
    ) { vector, angle -> LieTwist2d(vector, angle) }
    val poseArb: Arb<LiePose2d> = Arb.bind(
        vectorArb,
        rotationArb
    ) { vector, heading -> LiePose2d(vector, heading) }

    "rotation composed with the identity should equal itself" {
        forAll(rotationArb) { rot ->
            rot * LieRotation2d() == rot
        }
    }

    "rotation composed with its inverse should equal the identity rotation" {
        forAll(rotationArb) { rotation ->
            val inverse = rotation.inverse()
            rotation * inverse == LieRotation2d()
        }
    }

    "rotation composition should be associative" {
        forAll(rotationArb, rotationArb, rotationArb) { rot1, rot2, rot3 ->
            (rot1 * rot2) * rot3 == rot1 * (rot2 * rot3)
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
            pose * LiePose2d() == pose
        }
    }

    "pose composed with its inverse should equal the identity pose" {
        forAll(poseArb) { pose ->
            val inverse = pose.inverse()
            pose * inverse == LiePose2d()
        }
    }

    "pose composition should be associative" {
        forAll(poseArb, poseArb, poseArb) { pose1, pose2, pose3 ->
            (pose1 * pose2) * pose3 == pose1 * (pose2 * pose3)
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

    "pose exp and log should invert each other" {
        forAll(twistArb) { twist ->
            LiePose2d.exp(twist).log() == twist
        }
    }
})