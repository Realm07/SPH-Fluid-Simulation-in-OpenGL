#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <random>
#include <chrono>
#include <iomanip>
#include <imgui.h>
#include <omp.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#define GL_SILENCE_DEPRECATION
#include <cmath>

#ifndef M_PI
#define M_PI 3.14
#endif

constexpr auto SCREEN_WIDTH = 1280;
constexpr auto SCREEN_HEIGHT = 720; 

int rows = 20;
int cols = rows;
int radius = 3;
float spacingX = 2 * radius + 1;
float spacingY = 2 * radius + 1;
float width = SCREEN_WIDTH / 2 - (((cols + 1) * spacingX) / 2);
float height = SCREEN_HEIGHT / 2 - (((cols + 1) * spacingX) / 2);
static float smoothingRadius = 50.0f;
float targetDensity = 0.0001f;
float pressureMultiplier = 0.01f;
static float gravity = 0.000f;
float dampingFactor = 0.6;
float multiplicativeFactor = 1.0f;
bool show_directional_lines = false;
int nbFrames = 0;
double lastTime = glfwGetTime();
int iterations = 0;
float boundsSizeX = SCREEN_WIDTH;
float boundsSizeY = SCREEN_HEIGHT;

float fast_hypot(float x, float y) {
    const float a = std::abs(x);
    const float b = std::abs(y);
    return (a > b) ? (a + 0.960f * b) : (b + 0.960f * a);
}

class Vector2 {
public:
    float X;
    float Y;

    Vector2(float x, float y) : X(x), Y(y) {}

    float magnitude() const {
        return fast_hypot(X, Y);
    }
    Vector2 operator-(const Vector2& other) const {
        return Vector2(X - other.X, Y - other.Y);
    }
    Vector2 operator/(const float& scalar) const {
        return Vector2(X / scalar, Y / scalar);
    }
    Vector2 operator*(const float& right) const {
        return Vector2(X * right, Y * right);
    }
    Vector2 operator-(const float& scalar) const {
        return Vector2(X - scalar, Y - scalar);
    }
    Vector2 operator*(const Vector2& other) const {
        return Vector2(X * other.X, Y * other.Y);
    }
    Vector2 operator-() const {
        return Vector2(-X, -Y);
    }
    Vector2 operator+(const Vector2& other) const {
        return Vector2(X + other.X, Y + other.Y);
    }
    Vector2& operator+=(const Vector2& vec) {
        X += vec.X;
        Y += vec.Y;
        return *this;
    }
    static Vector2 Zero() {
        return Vector2(0.0f, 0.0f);
    }
};
void Update()
{
	velocity += Vector2.down * gravity;
	position += velocity;
	ResolveCollisions();
	ball.draw;
}

class Ball {
public:
    float x;
    float y;
    float vx;
    float vy;
    float radius;

    Ball(float x, float y, float vx, float vy, float radius) : x(x), y(y), vx(vx), vy(vy), radius(radius) {}


    void resolveCollision(float boundsSizeX, float boundsSizeY) {
        if (x + radius > boundsSizeX) {
            x = boundsSizeX - radius;
            vx *= -1 * dampingFactor;
        }
        else if (x - radius < 0) {
            x = radius;
            vx *= -1 * dampingFactor;
        }

        if (y + radius > boundsSizeY) {
            y = boundsSizeY - radius;
            vy *= -1 * dampingFactor;
        }
        else if (y - radius < 0) {
            y = radius;
            vy *= -1 * dampingFactor;
        }
    }

    void draw() const {
        drawCircle(x, y, 0, radius, 64);
    }

public:
    void drawCircle(GLfloat x, GLfloat y, GLfloat z, GLfloat radius, GLint numberOfSides) const {
        GLfloat twicePi = 2.0f * M_PI;

        GLfloat* allCircleVertices = new GLfloat[(numberOfSides + 2) * 3];
        allCircleVertices[0] = x;
        allCircleVertices[1] = y;
        allCircleVertices[2] = z;

        for (int i = 1; i < numberOfSides + 2; i++)
        {
            allCircleVertices[i * 3] = x + (radius * cos(i * twicePi / numberOfSides));
            allCircleVertices[i * 3 + 1] = y + (radius * sin(i * twicePi / numberOfSides));
            allCircleVertices[i * 3 + 2] = z;
        }

        glColor3f(0.0f, 0.8f, 1.0f);
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_FLOAT, 0, allCircleVertices);
        glDrawArrays(GL_TRIANGLE_FAN, 0, numberOfSides + 2);
        glDisableClientState(GL_VERTEX_ARRAY);

        delete[] allCircleVertices;
    }
};

static void drawCircle(float x, float y, float radius, int numSegments) {
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < numSegments; i++) {
        float theta = 2.0f * M_PI * static_cast<float>(i) / static_cast<float>(numSegments);
        float dx = radius * cosf(theta);
        float dy = radius * sinf(theta);
        glVertex2f(x + dx, y + dy);
    }
    glEnd();
}

std::vector<Ball> balls;

int main(void)
{
    std::cout << "Initializing GLFW" << std::endl;
    GLFWwindow* window;

    if (!glfwInit())
    {
        return -1;
    }
    std::cout << "Initializing Window with width " << SCREEN_WIDTH << " and height " << SCREEN_HEIGHT << std::endl;
    window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Simulation", NULL, NULL);

    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    std::cout << "Making window current Context" << std::endl;
    glfwMakeContextCurrent(window);
    std::cout << "Setting up viewport and other initializations" << std::endl;
    glViewport(0.0f, 0.0f, SCREEN_WIDTH, SCREEN_HEIGHT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0, 1);
    glMatrixMode(GL_MODELVIEW);

    //IMGUI
    std::cout << "Initializing IMGUI" << std::endl;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    const char* glsl_version = "#version 130";
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;

    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);
    bool show_smoothing_radius = false;

    glfwSwapInterval(1);

    std::mt19937 mt{ std::random_device{}() };
    std::uniform_real_distribution<float> distX(0.0f, SCREEN_WIDTH);
    std::uniform_real_distribution<float> distY(0.0f, SCREEN_HEIGHT);


    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j) {
            float x = distX(mt);
            float y = distY(mt);
            balls.emplace_back(x, y, 0.0, 0.0, radius);
        }
    }

    //Grid

    //for (int i = 0; i < rows; ++i) 
    //{
    //    for (int j = 0; j < cols; ++j) {
    //        float x = (j + 1) * spacingX;
    //        float y = (i + 1) * spacingY;
    //        //spacingX += 0.3;
    //        //spacingY += 0.3;
    //        balls.emplace_back(width + x,height + y, 0.0, 0.0, radius);
    //    }
    //}
    std::cout << "Entering window loop" << std::endl;

    while (!glfwWindowShouldClose(window))
    {
        glClear(GL_COLOR_BUFFER_BIT);

        std::vector<Vector2> positions;
        std::vector<Vector2> velocity;
        positions.reserve(balls.size()); // Reserve memory for expected number of elements

        #pragma omp parallel for
        for (auto& ball : balls) {
            ball.draw();
            ball.applyGravity(gravity);
            positions.push_back(Vector2(ball.x, ball.y));
        }

#pragma omp parallel for
        /*for (int i = 0; i < positions.size(); ++i) {
            densities[i] = calculator.CalculateDensity(positions[i]);
        }*/

        for (auto& ball : balls) {
            ball.updatePosition();
            ball.resolveCollision(boundsSizeX, boundsSizeY);
        }

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        {
            static int ballRadius = radius;
            ImGui::Begin("Config");
            ImGui::SliderInt("Ball Radius", &ballRadius, 1, 70);
            ImGui::End();

            if (radius != ballRadius) {
                radius = ballRadius;
                for (auto& ball : balls) {
                    ball.radius = radius;

                }
            }
        }
        /*if (show_smoothing_radius)
        {
            for (const auto& pos : positions) {
                drawBounds(pos, smoothingRadius);
            }

        }
        if (show_directional_lines)
        {
            show_directional_lines = true;

        }*/
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}