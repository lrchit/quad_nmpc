#include "orientation_tools.h"

Quat ori::quatProduct(Quat& q1, Quat& q2) {
    double r1 = q1[0];
    double r2 = q2[0];
    Vec3 v1(q1[1], q1[2], q1[3]);
    Vec3 v2(q2[1], q2[2], q2[3]);

    double r = r1 * r2 - v1.dot(v2);
    Vec3 v = r1 * v2 + r2 * v1 + v1.cross(v2);

    Quat q(r, v[0], v[1], v[2]);
    return q;
}

void ori::quaternionToso3(const Quat quat, Vec3& so3) {
    so3[0] = quat[1];
    so3[1] = quat[2];
    so3[2] = quat[3];

    double theta =
        2.0 * asin(sqrt(so3[0] * so3[0] + so3[1] * so3[1] + so3[2] * so3[2]));
    
    if (fabs(theta) < 0.0000001) {
        so3.setZero();
        return;
    }
    so3 /= sin(theta / 2.0);
    so3 *= theta; 
}

// rotation matrix is global to lacal
Mat3 ori::quaternionToRotationMatrix(const Quat& q) {
    double e0 = q(0);
    double e1 = q(1);
    double e2 = q(2);
    double e3 = q(3);

    Mat3 R;

    R << 1 - 2 * (e2 * e2 + e3 * e3), 2 * (e1 * e2 - e0 * e3),
        2 * (e1 * e3 + e0 * e2), 2 * (e1 * e2 + e0 * e3),
        1 - 2 * (e1 * e1 + e3 * e3), 2 * (e2 * e3 - e0 * e1),
        2 * (e1 * e3 - e0 * e2), 2 * (e2 * e3 + e0 * e1),
        1 - 2 * (e1 * e1 + e2 * e2);
    R.transposeInPlace();
    return R;
}

// be carful with rot order
Quat ori::so3ToQuat(Vec3& so3) {
    Quat quat;

    double theta = sqrt(so3[0] * so3[0] + so3[1] * so3[1] + so3[2] * so3[2]);

    if (fabs(theta) < 1.e-6) {
        quat.setZero();
        quat[0] = 1.;
        return quat;
    }
    quat[0] = cos(theta / 2.);
    quat[1] = so3[0] / theta * sin(theta / 2.);
    quat[2] = so3[1] / theta * sin(theta / 2.);
    quat[3] = so3[2] / theta * sin(theta / 2.);
    return quat;
}


/*!
 * Convert a quaternion to RPY.  Uses ZYX order (yaw-pitch-roll), but returns
 * angles in (roll, pitch, yaw).
 */
Vec3 ori::quatToRPY(const Quat& q) {
  Vec3 rpy;
  double as = std::min(-2. * (q[1] * q[3] - q[0] * q[2]), .99999);
  rpy(2) =
      std::atan2(2 * (q[1] * q[2] + q[0] * q[3]),
                 q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3]);
  rpy(1) = std::asin(as);
  rpy(0) =
      std::atan2(2 * (q[2] * q[3] + q[0] * q[1]),
                 q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3]);
  return rpy;
}

Quat ori::rpyToQuat(const Vec3& rpy) {

  Mat3 R = ori::rpyToRotMat(rpy);
  Quat q = ori::rotationMatrixToQuaternion(R);
  return q;
}

/*!
 * Compute rotation matrix for coordinate transformation. Note that
 * coordinateRotation(CoordinateAxis:X, .1) * v will rotate v by -.1 radians -
 * this transforms into a frame rotated by .1 radians!.
 */
Mat3 ori::coordinateRotation(ori::CoordinateAxis axis, double theta) {

  double s = std::sin(theta);
  double c = std::cos(theta);

  Mat3 R;

  if (axis == ori::CoordinateAxis::X) {
    R << 1, 0, 0, 0, c, s, 0, -s, c;
  } else if (axis == ori::CoordinateAxis::Y) {
    R << c, 0, -s, 0, 1, 0, s, 0, c;
  } else if (axis == ori::CoordinateAxis::Z) {
    R << c, s, 0, -s, c, 0, 0, 0, 1;
  }

  return R;
}

/*!
 * Go from rpy to rotation matrix.
 */
Mat3 ori::rpyToRotMat(const Vec3& v) {
  Mat3 m = ori::coordinateRotation(ori::CoordinateAxis::X, v[0]) *
                                    ori::coordinateRotation(ori::CoordinateAxis::Y, v[1]) *
                                    ori::coordinateRotation(ori::CoordinateAxis::Z, v[2]); //[todo] maybe error
  // Mat3 m = ori::coordinateRotation(ori::CoordinateAxis::Z, v[2]) *
  //                                   ori::coordinateRotation(ori::CoordinateAxis::Y, v[1]) *
  //                                   ori::coordinateRotation(ori::CoordinateAxis::X, v[0]);
  return m;
}

/*!
 * Convert a coordinate transformation matrix to an orientation quaternion.
 */
Quat ori::rotationMatrixToQuaternion(const Mat3& r1) {
  Quat q;
  Mat3 r = r1.transpose();
  double tr = r.trace();
  if (tr > 0.0) {
    double S = sqrt(tr + 1.0) * 2.0;
    q(0) = 0.25 * S;
    q(1) = (r(2, 1) - r(1, 2)) / S;
    q(2) = (r(0, 2) - r(2, 0)) / S;
    q(3) = (r(1, 0) - r(0, 1)) / S;
  } else if ((r(0, 0) > r(1, 1)) && (r(0, 0) > r(2, 2))) {
    double S = sqrt(1.0 + r(0, 0) - r(1, 1) - r(2, 2)) * 2.0;
    q(0) = (r(2, 1) - r(1, 2)) / S;
    q(1) = 0.25 * S;
    q(2) = (r(0, 1) + r(1, 0)) / S;
    q(3) = (r(0, 2) + r(2, 0)) / S;
  } else if (r(1, 1) > r(2, 2)) {
    double S = sqrt(1.0 + r(1, 1) - r(0, 0) - r(2, 2)) * 2.0;
    q(0) = (r(0, 2) - r(2, 0)) / S;
    q(1) = (r(0, 1) + r(1, 0)) / S;
    q(2) = 0.25 * S;
    q(3) = (r(1, 2) + r(2, 1)) / S;
  } else {
    double S = sqrt(1.0 + r(2, 2) - r(0, 0) - r(1, 1)) * 2.0;
    q(0) = (r(1, 0) - r(0, 1)) / S;
    q(1) = (r(0, 2) + r(2, 0)) / S;
    q(2) = (r(1, 2) + r(2, 1)) / S;
    q(3) = 0.25 * S;
  }
  return q;
}

