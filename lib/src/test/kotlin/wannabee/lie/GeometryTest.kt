/*
 * This Kotlin source file was generated by the Gradle 'init' task.
 */
package wannabee.lie

import io.kotest.core.spec.style.StringSpec
import io.kotest.property.Arb
import io.kotest.property.arbitrary.double
import io.kotest.property.forAll
import kotlin.math.PI

class GeometryTest: StringSpec({
    "pose composed with the identity should equal itself" {
        forAll(Arb.double(-1000.0..1000.0), Arb.double(-1000.0..1000.0), Arb.double(-1000.0..1000.0)) { x, y, heading ->
            val pose = LiePose2d(x, y, heading)
            pose * LiePose2d.identity() == pose
        }
    }

    "pose composed with its inverse should equal the identity pose" {
        forAll(Arb.double(-1000.0..1000.0), Arb.double(-1000.0..1000.0), Arb.double(-1000.0..1000.0)) { x, y, heading ->
            val pose = LiePose2d(x, y, heading)
            val inverse = pose.inverse()
            pose * inverse == LiePose2d.identity()
        }
    }
})