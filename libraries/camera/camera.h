/*

* Copyright (c) 2007, Carnegie Mellon University
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of Carnegie Mellon University, nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  Code author: Jernej Barbic, 2003-2007
  Research: Jernej Barbic, Doug L. James
  Version: 1.0

  A "spherical" OpenGL camera class, with the ability to set the OpenGL modelview transformation matrix.

  A "spherical" camera is a camera that is located at location (R, Phi, Theta), in spherical coordinates, away from some focus position. It is pointed towards the focus position. It is useful for orbiting a fixed location (focus position) in interactive applications.

  cameraPosition = focusPosition + 
   [ R * cos(Phi) * cos (Theta), R * sin(Theta), -R * sin(Phi) * cos (Theta) ]

*/

#ifndef _CAMERA_H_
#define _CAMERA_H_

#ifdef WIN32
  #include <windows.h>
#endif

class SphericalCamera
{
public:

  // R is the camera distance from the focus position 
  // Phi, Theta are camera longitude and lattitude (in radians)
  // focusPosition is the 3D camera focus position (e.g., [0, 0, 0])
  // up is the 3D up vector (e.g., [0, 1, 0])
  // movementSensitivity controls how much camera moves during the camera move commands
  // camera2WorldScalingFactor controls the scaling between the "camera" coordinate systems and "world" coordinate systems (only needed with the camera<-->world conversion routines)
  SphericalCamera(double R, double Phi, double Theta, double * focusPosition, double * up, double movementSensitivity=1.0, double camera2WorldScalingFactor=1.0);

  // commands to move the camera
  void MoveRight(double);
  void MoveUp(double);
  void ZoomIn(double); // negative value of parameter means zoom out
  void MoveIn(double);
  void MoveFocusRight(double); // move the camera focus, relative to the current focus
  void MoveFocusUp(double);

  void Look(); // calls gluLookAt to set the OpenGL modelview matrix, corresponding to the current camera position, focus, and up vector

  double GetRadius() { return R; }
  double GetPhi() { return Phi; }
  double GetTheta() { return Theta; }

  void SetPosition(double r, double phi, double theta);
  void SetFocusPosition(double focusPosition_[3]);
  void GetFocusPosition(double focusPosition_[3]);

  // get the camera position in the world coordinate system
  void GetAbsWorldPosition(double & x, double & y, double & z);

  void SetDefault(); // set current camera position as the default position
  void Reset(); // reset camera position to the default position
  
  // save/load current camera position to/from a file
  void SavePosition(const char * filename);
  void LoadPosition(const char * filename);

  // push/pop camera position onto an internal stack
  void PushPosition();
  void PopPosition();

  // determines good camera viewing parameters for an object centered at [centerX, centerY, centerZ], of radius "modelRadius"
  static void DetermineCameraParameters(double centerX, double centerY, double centerZ, double modelRadius, double * focusX, double * focusY, double * focusZ, double * cameraRadius, double * zNear, double * zFar);

  // get observer's ear positions, assuming given ear separation distance
  // observer is located at the camera and faces the focus position
  // (useful with audio simulations)
  void GetStereoPosition(double earSeparation,double * leftEarX, double * leftEarY, double * leftEarZ, double * rightEarX, double * rightEarY, double * rightEarZ);

  // === routines below this point are advanced ===

  double GetCamera2WorldScalingFactor() { return camera2WorldScalingFactor; }
  void SetCamera2WorldScalingFactor(double factor) { camera2WorldScalingFactor = factor; }

  // sets the world-coordinate position corresponding to zero position in the "camera" coordinate system
  inline void SetOrigin(double * origin_) { origin[0] = origin_[0]; origin[1] = origin_[1]; origin[2] = origin_[2]; }

  // transform a location specified in the camera coordinate system into the world coordinate system
  void CameraVector2WorldVector2D(double * c, double * w);
  void CameraVector2WorldVector2D(double cx, double cy, double cz, double * w);

  // same as above, except that it doesn't take into account the position of the camera
  // useful for transforming velocities
  void CameraVector2WorldVector_OrientationOnly2D(double * c, double * w);
  void CameraVector2WorldVector_OrientationOnly2D(double cx, double cy, double cz, double * w);
  void CameraVector2WorldVector_OrientationOnly3D(double cx, double cy, double cz, double * w); // uses a coordinate system where z-axis is normal to the camera sphere

  // converts a rotation in camera system to the corresponding rotation in the world coordinate system
  void CameraRotation2WorldRotation2D(float * c, float * w);
  void CameraRotation2WorldRotation2D(double * c, double * w);

  // converts 4x4 row-major camera-space transform matrix to 4x4 row-major world-space transform matrix
  void CameraTransform2WorldTransform2D(double * c, double * w); 
  void CameraTransform2WorldTransform2D(float * c, float * w); 
  void CameraTransform2WorldTransform2D_ColumnMajor(float * c, float * w); 
  void CameraTransform2WorldTransform2D_NoScaling_ColumnMajor(float * c, float * w);

  // transform vector given in the world coordinate system into the camera coordinate system
  // assumes the camera coordinate system and the world coordinate system origins coincide
  // useful for transforming forces and torques from world coordinate system back to the user(=camera) coordinate system
  void WorldVector2CameraVector_OrientationOnly2D(double * w, double * c);
  void WorldVector2CameraVector_OrientationOnly2D(float * w, float * c);
  void WorldVector2CameraVector_OrientationOnly2D(double w0, double w1, double w2, double * c);

  void Get2DAxes(double xAxis2D[2], double yAxis2D[2]);
  void Get3DAxes(double xAxis3D[3], double yAxis3D[3], double zAxis3D[3]);

  void ComputeCameraPosition(); // the user normally never needs to call this routine
  void ComputeLocalCoordinateSystem(); // the user normally never needs to call this routine

protected:

  // the default values
  double R0, Phi0, Theta0;

  double focusPosition[3]; // where the camera is aiming at, in world coordinate frame
  double up[3]; // the up vector, in world coordinate frame
  double cameraPosition[3]; // position of the camera, in world coordinate frame

  double R,Phi,Theta; // the current values

  // x_world = cameraMatrix * x_camera
  //Mat3f cameraMatrix;
  //Mat3f cameraMatrixT;

  double buf[9];

  double xAxis[3], yAxis[3], zAxis[3]; // the camera coordinate system axes; xAxis is parallel to xz-plane, zAxis points towards focus position

  double xAxis2D[3], yAxis2D[3]; // the simplified view axes (for 2D viewing transformation only)
     // these are actually 2d vectors, third components is only used for easy normalizing via NORMALIZE macro

  double movementSensitivity;
  double camera2WorldScalingFactor; // dx(world) = dx(camera) * camera2WorldScalingFactor

  double origin[3]; // only used to transform from camera to world coordinate system

  struct {
    double R,Phi,Theta;
    double focusPosition[3];
  } savedState;
};

#endif

