#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <math.h>
#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>
#include <iostream>
#include <sstream>
#ifndef M_PI
#define M_PI 3.1415
#endif

#define SCREEN_WIDTH 1280
#define SCREEN_HEIGHT 720

void drawCircle(GLfloat x, GLfloat y, GLfloat z, GLfloat radius, GLint numberOfSides);

class Vector2 {
public:
    float X;
    float Y;

    Vector2(float x, float y) : X(x), Y(y) {}
}; 

void drawLine(float x1, float y1, float x2, float y2);
void drawBounds(Vector2 size);

float radius = 20;
float gravity = 0.0010;
float dampingFactor = 0.88;
Vector2 boundsSize(SCREEN_WIDTH - 20, SCREEN_HEIGHT - 20);
Vector2 position(SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2);
Vector2 velocity(0.1, 0);


float Sign(float value) {
    return (value > 0) ? 1 : ((value < 0) ? -1 : 0);
}

float Abs(float value) {
    return (value < 0) ? -value : value;
}

void Update() {
    velocity.Y -= gravity;
    position.X += velocity.X;
    position.Y += velocity.Y;
}
void ResolveCollisions() {
   
    float halfBoundsWidth = boundsSize.X / 2;
    float halfBoundsHeight = boundsSize.Y / 2;

    float left = SCREEN_WIDTH / 2 - halfBoundsWidth;
    float right = SCREEN_WIDTH / 2 + halfBoundsWidth;
    float top = SCREEN_HEIGHT / 2 + halfBoundsHeight;
    float bottom = SCREEN_HEIGHT / 2 - halfBoundsHeight;


    if (position.X - radius < left) {
        position.X = left + radius;
        velocity.X *= -1 * dampingFactor;
    }
    else if (position.X + radius > right) {
        position.X = right - radius;
        velocity.X *= -1 * dampingFactor;
    }

    if (position.Y - radius < bottom) {
        position.Y = bottom + radius;
        velocity.Y *= -1 * dampingFactor;
    }
    else if (position.Y + radius > top) {
        position.Y = top - radius;
        velocity.Y *= -1 * dampingFactor;
    }
}


int main(void)
{
    GLFWwindow* window;

    // Initialize the library
    if (!glfwInit())
    {
        return -1;
    }
    window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Hello World", NULL, NULL); //Windowed mode

    if (!window)
    {
        glfwTerminate();
        return -1;
    }


    glfwMakeContextCurrent(window);

    glViewport(0.0f, 0.0f, SCREEN_WIDTH, SCREEN_HEIGHT); // specifies the part of the window to which OpenGL will draw (in pixels), convert from normalised to pixels
    glMatrixMode(GL_PROJECTION); // projection matrix defines the properties of the camera that views the objects in the world coordinate frame. Here you typically set the zoom factor, aspect ratio and the near and far clipping planes
    glLoadIdentity(); // replace the current matrix with the identity matrix and starts us a fresh because matrix transforms such as glOrpho and glRotate cumulate, basically puts us at (0, 0, 0)
    glOrtho(0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0, 1); // essentially set coordinate system
    glMatrixMode(GL_MODELVIEW); // (default matrix mode) modelview matrix defines how your objects are transformed (meaning translation, rotation and scaling) in your world
    glLoadIdentity(); // same as above comment

    double lastTime = glfwGetTime();
    int nbFrames = 0;

    
    while (!glfwWindowShouldClose(window)) {
        double currentTime = glfwGetTime();
        nbFrames++;
        if (currentTime - lastTime >= 1.0) {
            std::stringstream ss;
            ss << "FPS: " << nbFrames;

            glfwSetWindowTitle(window, ss.str().c_str());

            nbFrames = 0;
            lastTime += 1.0;
        }
        glClear(GL_COLOR_BUFFER_BIT);

        Update();
        ResolveCollisions();
       
        drawCircle(position.X, position.Y, 0, radius, 64);
        drawBounds(boundsSize);
        
        glfwSwapBuffers(window);

        glfwPollEvents();
    }


    glfwTerminate();

    return 0;
}

void drawCircle(GLfloat x, GLfloat y, GLfloat z, GLfloat radius, GLint numberOfSides)
{
    int numberOfVertices = numberOfSides + 2;

    GLfloat twicePi = 2.0f * M_PI;

    GLfloat* circleVerticesX = new GLfloat[numberOfVertices];
    GLfloat* circleVerticesY = new GLfloat[numberOfVertices];
    GLfloat* circleVerticesZ = new GLfloat[numberOfVertices];

    circleVerticesX[0] = x;
    circleVerticesY[0] = y;
    circleVerticesZ[0] = z;

    for (int i = 1; i < numberOfVertices; i++)
    {
        circleVerticesX[i] = x + (radius * cos(i * twicePi / numberOfSides));
        circleVerticesY[i] = y + (radius * sin(i * twicePi / numberOfSides));
        circleVerticesZ[i] = z;
    }

    GLfloat* allCircleVertices = new GLfloat[numberOfVertices * 3];

    for (int i = 0; i < numberOfVertices; i++)
    {
        allCircleVertices[i * 3] = circleVerticesX[i];
        allCircleVertices[(i * 3) + 1] = circleVerticesY[i];
        allCircleVertices[(i * 3) + 2] = circleVerticesZ[i];
    }

    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, allCircleVertices);
    glDrawArrays(GL_TRIANGLE_FAN, 0, numberOfVertices);
    glDisableClientState(GL_VERTEX_ARRAY);
    
}
void drawLine(float x1, float y1, float x2, float y2) {
    glBegin(GL_LINES);
    glVertex2f(x1, y1);
    glVertex2f(x2, y2);
    glEnd();
}
void drawBounds(Vector2 size) {
    float halfWidth = size.X / 2;
    float halfHeight = size.Y / 2;

    float left = SCREEN_WIDTH / 2 - halfWidth;
    float right = SCREEN_WIDTH / 2 + halfWidth;
    float top = SCREEN_HEIGHT / 2 + halfHeight;
    float bottom = SCREEN_HEIGHT / 2 - halfHeight;

    glColor3f(1.0f, 1.0f, 1.0f); // White color
    drawLine(left, top, right, top);
    drawLine(right, top, right, bottom);
    drawLine(right, bottom, left, bottom);
    drawLine(left, bottom, left, top);
}

