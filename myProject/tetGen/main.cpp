// tetGen.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"



GLint winWidth=800,winHeight=600;
float alpha=M_PI_2*0.5f,beta=M_PI_2 * 0.3f;
float r=20;  //这个变量用于与alpha, beta协同确定Lx, Ly, Lz（也就是相机位置）的值
float fovDeg = 60;  //某角度你懂的
float fov = fovDeg / 180.0f * M_PI;  //某角度的rad制表示
float camSrc[3] = {0};
float camTar[3] = {0};
/**********鼠标控制需要使用的变量，目前感觉它并没有实际生效*******************/
float lastMouseX,lastMouseY;  
bool isLeftButtonDown;
bool isMiddleButtonDown;
bool isRightButtonDown;
/**********鼠标控制变量End*********************/

/*************键盘控制变量******************/
bool isFirstMove;
bool isAltKeyDown;

//用于控制是否显示左下角小方框
bool isShowScene = true;
//用于控制是否显示左下角方框中的坐标轴
bool isShowGrid = false;
/*************键盘控制变量End***************/

bool isComputing=false;
bool isEyeChanged = false;
bool isInitiate = true;

/*************** fem ******************/
using namespace RigFEM;
RigFEM::FEMSystem g_fem;

void init()
{
	glClearColor(0.5,0.5,0.5,0.5);
	glEnableClientState(GL_VERTEX_ARRAY);

	//initialize camera
	//相机的位置，由r, beta, alpha三个值共同确定。
	camSrc[0]=r*cos(beta)*cos(alpha);
	camSrc[1]=r*cos(beta)*sin(alpha);
	camSrc[2]=r*sin(beta);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(camSrc[0],camSrc[1],camSrc[2],camTar[0],camTar[1],camTar[2],0,0,1);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fovDeg,float(winWidth)/float(winHeight), 0.1f,1000.0f);

	// init data
	g_fem.init();
}

void displayFcn()
{
	glClearColor(1,1,1,1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (isShowGrid)
	{
		int nGrid = 10;
		glColor3f(0.7f,0.7f,0.7f);
		glBegin(GL_LINES);
		for (int ithGrid = -nGrid; ithGrid <= nGrid; ++ithGrid)
		{
			glVertex3f(-10, ithGrid, 0);
			glVertex3f(10, ithGrid, 0);
			glVertex3f(ithGrid, -10,0);
			glVertex3f(ithGrid, 10, 0);
		}

		for (int i = 0; i < 3; i++)
		{
			float clr[3] = {0};
			clr[i] = 1;
			glColor3fv(clr);
			glVertex3f(0,0,0);
			glVertex3fv(clr);
		}
		glEnd();
	}

	g_fem.show();

	glutSwapBuffers();
	glFlush();
}
void mouseFcn(GLint button,GLint action,GLint xMouse,GLint yMouse)
{
	int mod=glutGetModifiers();
	if(mod==GLUT_ACTIVE_ALT)
		isAltKeyDown=true;
	else
		isAltKeyDown=false;
	if(button==GLUT_LEFT_BUTTON)
	{		
		if(action==GLUT_DOWN)
		{
			isLeftButtonDown=true;
			isFirstMove=true;
		}
		else if(action==GLUT_UP)
		{
			isLeftButtonDown=false;
			isFirstMove=true;
		}		
	}
	if(button==GLUT_MIDDLE_BUTTON)
	{
		if(action==GLUT_DOWN)
		{
			isMiddleButtonDown = true;
			isFirstMove = true;
		}
		else if(action==GLUT_UP)
		{
			isMiddleButtonDown = false;
			isFirstMove=true;
		}
	}
	if(button==GLUT_RIGHT_BUTTON)
	{
		if(action==GLUT_DOWN)
		{
			isRightButtonDown=true;
			isFirstMove=true;
		}
		else if(action==GLUT_UP)
		{
			isRightButtonDown=false;
			isFirstMove=true;
		}
	}
}

void mouseMoveFcn(GLint xMouse,GLint yMouse)
{
	float deltaX;
	float deltaY;
	if(isLeftButtonDown&&isAltKeyDown)
	{
		if(isFirstMove)
		{
			lastMouseX=xMouse;
			lastMouseY=yMouse;
			isFirstMove=false;
		}
		else
		{		
			deltaX=xMouse-lastMouseX;
			deltaY=yMouse-lastMouseY;
			lastMouseX=xMouse;
			lastMouseY=yMouse;
			alpha-=deltaX*0.004f;//deltaX;
			beta+=deltaY*0.004f;//deltaY;
			camSrc[0]=r*cos(beta)*cos(alpha) + camTar[0];
			camSrc[1]=r*cos(beta)*sin(alpha) + camTar[1];
			camSrc[2]=r*sin(beta) + camTar[2];
		}
		isEyeChanged = true;
	}
	if(isMiddleButtonDown && isAltKeyDown)
	{
		if(isFirstMove)
		{
			lastMouseX=xMouse;
			lastMouseY=yMouse;
			isFirstMove=false;
		}
		else
		{
			deltaX=xMouse-lastMouseX;
			deltaY=yMouse-lastMouseY;
			lastMouseX=xMouse;
			lastMouseY=yMouse;

			float sinA = sin(alpha);
			float cosA = cos(alpha);
			float sinB = sin(beta);
			float cosB = cos(beta);

			float xDir[3] = {-sinA, cosA, 0};
			float yDir[3] = {-sinB*cosA, -sinB*sinA, cosB};

			for (int i = 0; i < 3; ++i)
			{
				float pMove = 0.01f * (-xDir[i] * deltaX + yDir[i] * deltaY);
				camTar[i] += pMove;
				camSrc[i] += pMove;
			}
		}
		isEyeChanged = true;
	}
	if(isRightButtonDown&&isAltKeyDown)
	{
		if(isFirstMove)
		{
			lastMouseX=xMouse;
			lastMouseY=yMouse;
			isFirstMove=false;
		}
		else
		{
			deltaX=xMouse-lastMouseX;
			lastMouseX=xMouse;
			r-=deltaX*0.05f;
			camSrc[0]=r*cos(beta)*cos(alpha) + camTar[0];
			camSrc[1]=r*cos(beta)*sin(alpha) + camTar[1];
			camSrc[2]=r*sin(beta) + camTar[2];
		}
		isEyeChanged = true;

	}

	if (isEyeChanged)
	{
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt(camSrc[0],camSrc[1],camSrc[2],camTar[0],camTar[1],camTar[2],0,0,1);
		glutPostRedisplay();
	}
}

void keyFcn(GLubyte key,GLint xMouse,GLint yMouse)
{
	switch (key)
	{
	case 'r':
		glutPostRedisplay();
		break;
	case 'g':
		isShowGrid = !isShowGrid;
		break;
	case 'e':
		isComputing = false;
		break;
	case 'v':
		isShowScene = !isShowScene; break;
	case 's':
		if (isComputing)
			g_fem.saveResult("paramRes.m");
		else
			g_fem.clearResult();
		isComputing = !isComputing;//g_fem.computeNextStep();
		break;
	case 't':
		g_fem.m_mesh.testValue();	
		//MathUtilities::testMath();
		break;
	case '[':
	case ']':
		TransformRig& rig =	g_fem.m_mesh.m_transRig;
		double dir = key == ']' ? 1.f : -1.f;
		Vec3d trans = rig.getTranslation() + Vec3d(0.0) * dir * 0;
		Vec3d scale = rig.getScale() + Vec3d(0.01) * dir * 0;
		Vec3d rotate= rig.getRotation() + Vec3d(0,0,1) * dir * 0.1;
		rig.setAllParam(trans, rotate, scale);
		g_fem.m_mesh.computeRig();
		break;
	}

	glutPostRedisplay();

}
void reshapeFcn(GLint Width,GLint Height)
{
	winWidth=Width;
	winHeight=Height;
	glViewport(0,0,Width, Height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fovDeg,float(winWidth)/float(winHeight), 0.1f,1000.0f);

}


void timerFcn(int v)
{
	glutTimerFunc(500, timerFcn, -1);
	glutPostRedisplay();
}

void idleFcn()
{
	if (isComputing)
	{
		g_fem.step();
	}
}
int main(int argc, char** argv)
{
	glutInit(&argc,(char**)argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(winWidth,winHeight);
	glutCreateWindow("WRay");


	init();

	glutIdleFunc(idleFcn);
	glutDisplayFunc(displayFcn);
	glutReshapeFunc(reshapeFcn);
	glutKeyboardFunc(keyFcn);
	glutMouseFunc(mouseFcn);
	glutMotionFunc(mouseMoveFcn);
	glutTimerFunc(500, timerFcn, -1);

	glutMainLoop();
	return 0;
}
