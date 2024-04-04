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
#define M_PI 3.15
#endif

constexpr auto SCREEN_WIDTH = 1280;
constexpr auto SCREEN_HEIGHT = 720;

int rows = 15;
int cols = rows;
int radius = 4;
int spacingX = 11;
int spacingY = spacingX;
float width = SCREEN_WIDTH / 2.0 - (((cols + 1.0) * spacingX) / 2.0);
float height = SCREEN_HEIGHT / 2.0 - (((cols + 1.0) * spacingX) / 2.0);
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

std::mt19937 mt{ std::random_device{}() };
std::uniform_real_distribution<float> distX(0.0f, SCREEN_WIDTH);
std::uniform_real_distribution<float> distY(0.0f, SCREEN_HEIGHT);

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

static float Sign(float value) {
    return (value > 0) ? 1 : ((value < 0) ? -1 : 0);
}
static float Abs(float value) {
    return (value < 0) ? -value : value;
}

static float SmoothingKernel(float radius, float dst)
{
    if (dst < radius)
    {
        float volume = (M_PI * pow(radius, 4)) / 6;
        float result = ((radius - dst) * (radius - dst) / volume);
        //clamp
        return (result >= 0.0f) ? result : 0.0f;

    }
    else return 0;

}
static float SmoothingKernelDerivative(float dst, float radius)
{
    if (dst < radius)
    {
        float scale = 12 / (pow(radius, 4) * M_PI);
        return (dst - radius) * scale;
    }
    else return 0;
}
static void glfw_error_callback(int error, const char* description)
{
    fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}

class Ball {
public:
    float x;
    float y;
    float vx;
    float vy;
    float radius;
    
    Ball(float x, float y, float vx, float vy, float radius) : x(x), y(y), vx(vx), vy(vy), radius(radius) {}

    void applyGravity(float gravity) {
        vy -= gravity;
    }

    void resolveCollision(float boundsSizeX, float boundsSizeY) {
        float minX = radius;
        float maxX = boundsSizeX - radius;
        float minY = radius;
        float maxY = boundsSizeY - radius;

        if (x < minX) {
            x = minX;
            vx *= -dampingFactor;
        }
        else if (x > maxX) {
            x = maxX;
            vx *= -dampingFactor;
        }

        if (y < minY) {
            y = minY;
            vy *= -dampingFactor;
        }
        else if (y > maxY) {
            y = maxY;
            vy *= -dampingFactor;
        }
    }

    void updatePosition() {
        x += vx;
        y += vy;
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

static void drawRadius(float x, float y, float z, float radius, int numSegments) {
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < numSegments; i++) {
        float theta = 2.0f * M_PI * static_cast<float>(i) / static_cast<float>(numSegments);
        float dx = radius * cosf(theta);
        float dy = radius * sinf(theta);
        glVertex3f(x + dx, y + dy, 0);
    }
    glEnd();
}

void drawBounds(Vector2 position, float radius) {
    glColor3f(1.0f, 0.0f, 0.0f);
    //drawRadius(position.X, position.Y, 0, radius, 30);
    drawRadius(SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2, 0, radius, 30);
}

float CalculateDensity(const std::vector<Vector2>& positions, float smoothingRadius, int index) {
    float density = 0;
    const float mass = 1;
    Vector2 samplePoint(SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2);

    Vector2 position = positions[index];
    for (const auto& p : positions) {
        Vector2 offset = p - position;
        float dst = offset.magnitude();

        if (dst > smoothingRadius) continue;

        float influence = SmoothingKernel(smoothingRadius, dst);
        density += mass * influence;
    }

    return density;
}

void PreCalculateDensities(std::vector<float>& densities, const std::vector<Vector2>& positions, float smoothingRadius) {
    if (densities.size() != positions.size()) {
        densities.resize(positions.size());
    }
    for (size_t i = 0; i < positions.size(); ++i) {
        densities[i] = CalculateDensity(positions, smoothingRadius, i);
    }
}

static float ConvertDensityToPressure(float density)
{
    float densityError = density - targetDensity;
    float pressure = densityError * pressureMultiplier;
    return pressure;
}

float CalculateSharedPressure(float densityA, float densityB)
{
    float pressureA = ConvertDensityToPressure(densityA);
    float pressureB = ConvertDensityToPressure(densityB);
    return (pressureA + pressureB) / 2;
};

static Vector2 getRandomDir() {
    float angle = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 2.0f * M_PI;
    return Vector2(cos(angle), sin(angle));
}

Vector2 randomDir = getRandomDir();

static Vector2 CalculatePressureForce(int particleIndex, const std::vector<Vector2>& positions, const std::vector<float>& densities, float smoothingRadius)
{
    Vector2 pressureForce = Vector2::Zero();
    Vector2 particlePosition = positions[particleIndex];
    float particleDensity = densities[particleIndex];

#pragma omp parallel for reduction(+:pressureForce)
    for (int otherParticleIndex = 0; otherParticleIndex < positions.size(); otherParticleIndex++)
    {
        if (particleIndex == otherParticleIndex) continue;

        Vector2 offset = positions[otherParticleIndex] - particlePosition;
        float dst = offset.magnitude();

        if (dst > smoothingRadius) continue;

        Vector2 dir = (dst == 0) ? randomDir : offset / dst;
        float slope = SmoothingKernelDerivative(dst, smoothingRadius);
        float density = densities[otherParticleIndex];
        float sharedPressure = CalculateSharedPressure(density, particleDensity);
        Vector2 scaledDir = dir * sharedPressure;
        pressureForce += scaledDir * slope / density;

        if (show_directional_lines == true) 
        {
            glBegin(GL_LINES);
            glColor3f(1.0f, 0.0f, 1.0f); 
            glVertex2f(particlePosition.X, particlePosition.Y);
            glVertex2f(positions[otherParticleIndex].X, positions[otherParticleIndex].Y);
            glEnd();
        }
    }
    return pressureForce;
}

float CalculateDensity(const std::vector<Vector2>& positions, float smoothingRadius)
{
    float density = 0;
    const float mass = 1;
    Vector2 samplePoint(SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2);

    for (const auto& p : positions) {
        Vector2 offset = p - samplePoint;
        float dst = offset.magnitude();

        if (dst > smoothingRadius) continue;

        float influence = SmoothingKernel(smoothingRadius, dst);
        density += mass * influence;
    }

    return density;
}

void initializeBalls(float spacingX, float spacingY, int rows, int cols, int radius, float width, float height, std::vector<Ball>& balls) {
    balls.clear();
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            float x = (j + 1) * spacingX;
            float y = (i + 1) * spacingY;
            balls.emplace_back(width + x, height + y, 0.0, 0.0, radius);
        }
    }
}

void updateBallPositions(float spacingX, float spacingY, int rows, int cols, int radius, float width, float height, std::vector<Ball>& balls) {
    balls.clear();
    initializeBalls(spacingX, spacingY, rows, cols, radius, width, height, balls);
}

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
    std::vector<Ball> balls;
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j) {
            float x = distX(mt);
            float y = distY(mt);
            balls.emplace_back(x, y, 0.0, 0.0, radius);
        }
    }

    /*for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            float x = (j + 1) * spacingX;
            float y = (i + 1) * spacingY;
            balls.emplace_back(width + x, height + y, 0.0, 0.0, radius);
        }
    }*/
    std::cout << "Entering window loop" << std::endl;
    
    
    while (!glfwWindowShouldClose(window))
    {
        glClear(GL_COLOR_BUFFER_BIT);

        std::vector<float> densities;
        std::vector<Vector2> positions;
        positions.reserve(balls.size());
      
        #pragma omp parallel for
        for (auto& ball : balls) {
           
            ball.applyGravity(gravity);
            positions.push_back(Vector2(ball.x, ball.y));
        }

        
        float density = CalculateDensity(positions, smoothingRadius);
        densities.resize(positions.size());

        PreCalculateDensities(densities, positions, smoothingRadius);

        #pragma omp parallel for
        for (int i = 0; i < balls.size(); ++i) {
            Vector2 pressureForce = CalculatePressureForce(i, positions, densities, smoothingRadius);
            Vector2 pressureAcceleration = pressureForce / densities[i];
            balls[i].vx += pressureAcceleration.X * multiplicativeFactor;
            balls[i].vy += pressureAcceleration.Y * multiplicativeFactor;
        }

        for (auto& ball : balls) {
            ball.updatePosition();
            ball.resolveCollision(boundsSizeX, boundsSizeY);
            ball.draw();
        }
        
        //width = SCREEN_WIDTH / 2 - (((cols + 1) * spacingX) / 2);
        //height = SCREEN_HEIGHT / 2 - (((cols + 1) * spacingY) / 2);
        //updateBallPositions(spacingX, spacingY, rows, cols, radius, width, height, balls);
        
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        {
            static int ballRadius = radius;
            ImGui::Begin("Config");
            ImGui::SliderFloat("Gravity", &gravity, 0.0000f, 0.10f);
            ImGui::SliderInt("Ball Radius", &ballRadius, 1, 70);
            ImGui::Text("Density: %.8f", density);
            ImGui::SliderFloat("Smoothing Radius", &smoothingRadius, 0.05f, 500.0f);
            ImGui::SliderInt("Spacing X", &spacingX, 0.0, 150.0);
            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
            ImGui::SliderFloat("Pressure Multiplier", &pressureMultiplier, 0.01f, 100.0f);
            //ImGui::SliderFloat("Multiplicative Factor", &multiplicativeFactor, 1.0f, 100.0f);
            //ImGui::SliderFloat("Target Density", &targetDensity, 0.0000001f, 1.0f);
            ImGui::Checkbox("Show Smoothing Radius", &show_smoothing_radius);
            ImGui::SameLine();
            ImGui::Checkbox("Show Directional Lines", &show_directional_lines);
            ImGui::End();

            if (radius != ballRadius) {
                radius = ballRadius;
                for (auto& ball : balls) {
                    ball.radius = radius;

                }
            }
            if (spacingY != spacingX)
            {
                spacingY = spacingX;
            }
        }
        if (show_smoothing_radius)
        {
            for (const auto& pos : positions) {
                drawBounds(pos, smoothingRadius);
            }

        }
        if (show_directional_lines)
        {
            show_directional_lines = true;
            
        }
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}

static void drawLine(float x1, float y1, float x2, float y2) {
    glBegin(GL_LINES);
    glVertex2f(x1, y1);
    glVertex2f(x2, y2);
    glEnd();
}
//void drawBounds(Vector2 size) {
//    float halfWidth = size.X / 2;
//    float halfHeight = size.Y / 2;
//
//    float left = SCREEN_WIDTH / 2 - halfWidth;
//    float right = SCREEN_WIDTH / 2 + halfWidth;
//    float top = SCREEN_HEIGHT / 2 + halfHeight;
//    float bottom = SCREEN_HEIGHT / 2 - halfHeight;
//
//    glColor3f(1.0f, 1.0f, 1.0f); // White color
//    drawLine(left, top, right, top);
//    drawLine(right, top, right, bottom);
//    drawLine(right, bottom, left, bottom);
//    drawLine(left, bottom, left, top);
//}
