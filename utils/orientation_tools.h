#ifndef ORIENTATION_TOOLS_H
#define ORIENTATION_TOOLS_H

#include <cmath>
#include "cppTypes.h"

namespace ori{
    enum class CoordinateAxis { X, Y, Z };

    Quat quatProduct(Quat& q1, Quat& q2);
    Mat3 quaternionToRotationMatrix(const Quat& q);
    Vec3 quatToRPY(const Quat& q);

    Quat so3ToQuat(Vec3& so3);
    void quaternionToso3(const Quat quat, Vec3& so3);

    Quat rpyToQuat(const Vec3& rpy);

    Mat3 coordinateRotation(ori::CoordinateAxis axis, double theta);
    Mat3 rpyToRotMat(const Vec3& v);
    Quat rotationMatrixToQuaternion(const Mat3& r1);

};
#endif