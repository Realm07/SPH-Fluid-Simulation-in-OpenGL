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
#include <cmath>
#include <tuple>
#include <algorithm>
#include <limits>
#include <utility>
#include <cstdint>
#include <array>
#include <unordered_set>
#include <set>
#include <execution>

constexpr float M_PI = 3.14;
constexpr GLfloat twicePi = 2.0f * M_PI;
int SCREEN_WIDTH = 1800;
int SCREEN_HEIGHT = 900;

int rows = 40;
int cols = rows;
int radius = 4;
constexpr int spacingX = 11;
constexpr int spacingY = spacingX;
float smoothingRadius = 66.0f;
float gravity = 0.000f;
constexpr float dampingFactor = 0.8f;
float targetDensity = 0.0003f;
float pressureMultiplier = 75.0f;
float stiffnessConstant = 1.0f;
float boundsSizeX = SCREEN_WIDTH;
float horizontalFactor = 0;
float verticalFactor = 0;
float boundsSizeY = SCREEN_HEIGHT;
float mouseRadius = 100;
constexpr GLint numberOfSides = 8;
float viscosityStrength = 15.0f;
float mouseStrength = 0.5f;
float sqrSmoothingRadius = smoothingRadius * smoothingRadius;
float lastFrame = 0.0f;
double mouseX, mouseY;
int lostParticles = 0;
float width = SCREEN_WIDTH / 2.0 - (((cols + 1.0) * spacingX) / 2.0);
float height = SCREEN_HEIGHT / 2.0 - (((cols + 1.0) * spacingX) / 2.0);

bool show_smoothing_radius = false;
bool show_directional_lines = false;
bool show_density_areas = false;
bool show_spatial_grid = false;
bool reset_to_random = true;
bool reset_to_grid = false;
bool show_mouse = true;
bool enableInteraction = false;

//Math
static float fast_hypot(float x, float y) {
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
    float sqrMagnitude() const {
        return X * X + Y * Y;
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
    Vector2& operator-=(const Vector2& other) {
        X -= other.X;
        Y -= other.Y;
        return *this;
    }
    static Vector2 Zero() {
        return Vector2(0.0f, 0.0f);
    }
    void normalize() {
        float mag = magnitude();
        if (mag != 0) {
            X /= mag;
            Y /= mag;
        }
    }
    static Vector2 down() {
        return Vector2(0.0f, -1.0f);
    }
};
static float Sign(float value) {
    return (value > 0) ? 1 : ((value < 0) ? -1 : 0);
}
static float Abs(float value) {
    return (value < 0) ? -value : value;
}

//Curves
static inline float SmoothingKernel(float radius, float dst)
{
    if (dst < radius)
    {
        float radius4 = radius * radius * radius * radius;
        float volume = (M_PI * radius4) / 6;
        float diff = radius - dst;
        float result = (diff * diff) / volume;

        return (result >= 0.0f) ? result : 0.0f;
    }
    else return 0;
}
static inline float SmoothingKernelDerivative(float dst, float radius)
{
    if (dst < radius)
    {
        float radius4 = radius * radius * radius * radius;
        float scale = 12 / (radius4 * M_PI);
        return (dst - radius) * scale;
    }
    else return 0;
}
static inline float ViscositySmoothingKernel(float dst, float radius)
{
    float radius8 = radius * radius * radius * radius * radius * radius * radius * radius;
    float volume = M_PI * radius8 / 4;
    float diff2 = radius * radius - dst * dst;
    return diff2 * diff2 * diff2 / volume;
}

//Randomization
static inline Vector2 getRandomDir() {
    float angle = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 2.0f * M_PI;
    return Vector2(cos(angle), sin(angle));
}
Vector2 randomDir = getRandomDir();

std::mt19937 mt{ std::random_device{}() };
std::uniform_real_distribution<float> distX(0.0f, SCREEN_WIDTH);
std::uniform_real_distribution<float> distY(0.0f, SCREEN_HEIGHT);

//GLFW
static void glfw_error_callback(int error, const char* description)
{
    fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}

//Particles
std::vector<Vector2> velocities;
class Ball {
public:
    Vector2 position;
    Vector2 velocity;
    float radius;

    Ball(const Vector2& position, const Vector2& velocity, float radius)
        : position(position), velocity(velocity), radius(radius) {}

    inline void draw() const {
        // Calculate the magnitude of the velocity
        float magnitude = sqrt(velocity.X * velocity.X + velocity.Y * velocity.Y);

        // Normalize based on expected max velocity
        float maxVelocity = 3.0f;
        float normalizedMagnitude = std::min(magnitude / maxVelocity, 1.0f);

        // Calculate the color based on the gradient
        float r, g, b;

        // Mapping normalizedMagnitude to gradient stops
        if (normalizedMagnitude <= 0.34f) {
            // Between rgba(12,0,202,1) and rgba(9,78,121,1)
            float t = normalizedMagnitude / 0.34f; // Transition factor
            r = (1 - t) * 12 + t * 9;
            g = (1 - t) * 0 + t * 78;
            b = (1 - t) * 202 + t * 121;
        }
        else if (normalizedMagnitude <= 0.52f) {
            // Between rgba(9,78,121,1) and rgba(38,169,57,1)
            float t = (normalizedMagnitude - 0.34f) / (0.52f - 0.34f); // Transition factor
            r = (1 - t) * 9 + t * 38;
            g = (1 - t) * 78 + t * 169;
            b = (1 - t) * 121 + t * 57;
        }
        else if (normalizedMagnitude <= 0.71f) {
            // Between rgba(38,169,57,1) and rgba(197,181,5,1)
            float t = (normalizedMagnitude - 0.52f) / (0.71f - 0.52f); // Transition factor
            r = (1 - t) * 38 + t * 197;
            g = (1 - t) * 169 + t * 181;
            b = (1 - t) * 57 + t * 5;
        }
        else {
            // Between rgba(197,181,5,1) and rgba(255,0,0,1)
            float t = (normalizedMagnitude - 0.71f) / (1.0f - 0.71f); // Transition factor
            r = (1 - t) * 197 + t * 255;
            g = (1 - t) * 181 + t * 0;
            b = (1 - t) * 5;
        }

        // Normalize the RGB values to [0, 1] for OpenGL
        r /= 255.0f;
        g /= 255.0f;
        b /= 255.0f;

        // Draw the ball with the calculated color
        drawCircle(position.X, position.Y, 0, radius, r, g, b);
    }


public:
    inline void drawCircle(GLfloat x, GLfloat y, GLfloat z, GLfloat radius, GLfloat r, GLfloat g, GLfloat b) const {
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

        glColor3f(r, g, b);
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_FLOAT, 0, allCircleVertices);
        glDrawArrays(GL_TRIANGLE_FAN, 0, numberOfSides + 2);
        glDisableClientState(GL_VERTEX_ARRAY);

        delete[] allCircleVertices;
    }
};

Vector2 obstacleCenter(1400.0f, 350.0f);
float obstacleWidth = 170.f;
float obstacleHeight = 375.0f;
// Draw the obstacle bounds

float minX = radius;
float maxX = boundsSizeX - radius;
float minY = radius;
float maxY = boundsSizeY - radius;

static inline void resolveCollisions(std::vector<Vector2>& positions, std::vector<Vector2>& velocities, float radius, float dampingFactor, float boundsSizeX, float boundsSizeY, const Vector2& obstacleCenter, float obstacleWidth, float obstacleHeight) {
    // Draw the main bounding box
    //glBegin(GL_LINE_LOOP);
    //glColor3f(1.0f, 1.0f, 1.0f);  // Set color to white
    //glVertex2f(0 + horizontalFactor, 0 + verticalFactor);
    //glVertex2f(boundsSizeX - horizontalFactor, 0 + verticalFactor);
    //glVertex2f(boundsSizeX - horizontalFactor, boundsSizeY - verticalFactor);
    //glVertex2f(0 + horizontalFactor, boundsSizeY - verticalFactor);
    //glEnd();

    float obstacleMinX = obstacleCenter.X - obstacleWidth / 2.0f;
    float obstacleMaxX = obstacleCenter.X + obstacleWidth / 2.0f;
    float obstacleMinY = obstacleCenter.Y - obstacleHeight / 2.0f;
    float obstacleMaxY = obstacleCenter.Y + obstacleHeight / 2.0f;

    glBegin(GL_LINE_LOOP);
    glColor3f(1.0f, 1.0f, 1.0f);
    glVertex2f(obstacleMinX, obstacleMinY);
    glVertex2f(obstacleMaxX, obstacleMinY);
    glVertex2f(obstacleMaxX, obstacleMaxY);
    glVertex2f(obstacleMinX, obstacleMaxY);
    glEnd();

    for (int i = 0; i < positions.size(); i++) {
        // Check collision with the bounding box
        if (positions[i].X < minX + horizontalFactor) {
            positions[i].X = minX + horizontalFactor;
            velocities[i].X *= -dampingFactor;
        }
        else if (positions[i].X > maxX - horizontalFactor) {
            positions[i].X = maxX - horizontalFactor;
            velocities[i].X *= -dampingFactor;
        }

        if (positions[i].Y < minY + verticalFactor) {
            positions[i].Y = minY + verticalFactor;
            velocities[i].Y *= -dampingFactor;
        }
        else if (positions[i].Y > maxY - verticalFactor) {
            positions[i].Y = maxY - verticalFactor;
            velocities[i].Y *= -dampingFactor;
        }

        // Check collision with the obstacle
        if (positions[i].X - radius < obstacleMaxX && positions[i].X + radius > obstacleMinX &&
            positions[i].Y - radius < obstacleMaxY && positions[i].Y + radius > obstacleMinY) {
            // Determine which side of the obstacle the particle is closest to and adjust accordingly
            float distToLeft = std::abs(positions[i].X - radius - obstacleMinX);
            float distToRight = std::abs(positions[i].X + radius - obstacleMaxX);
            float distToTop = std::abs(positions[i].Y + radius - obstacleMaxY);
            float distToBottom = std::abs(positions[i].Y - radius - obstacleMinY);

            // Resolve the collision with the nearest boundary of the obstacle
            if (distToLeft < distToRight && distToLeft < distToTop && distToLeft < distToBottom) {
                positions[i].X = obstacleMinX - radius;
                velocities[i].X *= -dampingFactor;
            }
            else if (distToRight < distToLeft && distToRight < distToTop && distToRight < distToBottom) {
                positions[i].X = obstacleMaxX + radius;
                velocities[i].X *= -dampingFactor;
            }
            else if (distToTop < distToLeft && distToTop < distToRight && distToTop < distToBottom) {
                positions[i].Y = obstacleMaxY + radius;
                velocities[i].Y *= -dampingFactor;
            }
            else if (distToBottom < distToLeft && distToBottom < distToRight && distToBottom < distToTop) {
                positions[i].Y = obstacleMinY - radius;
                velocities[i].Y *= -dampingFactor;
            }
        }

        
    }
}

//Drawing
static inline void drawRadius(float x, float y, float z, float radius, int numSegments) {
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < numSegments; i++) {
        float theta = 2.0f * M_PI * static_cast<float>(i) / static_cast<float>(numSegments);
        float dx = radius * cosf(theta);
        float dy = radius * sinf(theta);
        glVertex3f(x + dx, y + dy, 0);
    }
    glEnd();
}
static void drawBounds(float radius) {
    glColor3f(1.0f, 0.0f, 0.0f);
    //drawRadius(position.X, position.Y, 0, radius, 30);
    drawRadius(SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2, 0, radius, 30);
}
static inline void drawOutline(Vector2 position, float smoothingRadius) {
    glColor3f(0.0f, 1.0f, 0.0f);
    //drawRadius(position.X, position.Y, 0, radius, 30);
    drawRadius(position.X, position.Y, 0, smoothingRadius, 128);
}
static void drawLine(float x1, float y1, float x2, float y2) {
    glBegin(GL_LINES);
    glVertex2f(x1, y1);
    glVertex2f(x2, y2);
    glEnd();
}

//Density
static inline float CalculateDensity(const std::vector<Vector2>& positions, float smoothingRadius)
{
    float density = 0;
    const float mass = 1;
    Vector2 samplePoint(SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2);

    for (const auto& p : positions) {
        Vector2 offset = p - samplePoint;
        float sqrDst = offset.sqrMagnitude();

        if (sqrDst > sqrSmoothingRadius) continue;

        float influence = SmoothingKernel(smoothingRadius, std::sqrt(sqrDst));
        density += mass * influence;
    }

    return density;
}
static inline float CalculateDensity(int index, const std::vector<int>& particleIndices, const std::vector<Vector2>& positions, float smoothingRadius) {
    float density = 0;
    const float mass = 1;
    Vector2 position = positions[index];
    for (int pIndex : particleIndices) {
        Vector2 offset = positions[pIndex] - position;
        float sqrDst = offset.sqrMagnitude();

        if (sqrDst > sqrSmoothingRadius) continue;

        float influence = SmoothingKernel(smoothingRadius, std::sqrt(sqrDst));
        density += mass * influence;
    }

    return density;
}
//static void PreCalculateDensities(std::vector<float>& densities, const std::vector<int>& particleIndices, const std::vector<Vector2>& positions, float smoothingRadius) {
//    if (densities.size() != positions.size()) {
//        densities.resize(positions.size());
//    }
//    for (size_t i = 0; i < positions.size(); ++i) {
//        densities[i] = CalculateDensity(i, particleIndices, positions, smoothingRadius);
//    }
//}

//Forces
static inline float ConvertDensityToPressure(float density)
{
    return (density - targetDensity) * pressureMultiplier;
}
static inline float CalculateSharedPressure(float densityA, float densityB)
{
    float pressureA = ConvertDensityToPressure(densityA);
    float pressureB = ConvertDensityToPressure(densityB);
    return (pressureA + pressureB) / 2;
};
static inline Vector2 CalculatePressureForce(int index, const std::vector<int>& particleIndices, const std::vector<Vector2>& positions, const std::vector<float>& densities, float smoothingRadius)
{
    Vector2 pressureForce = Vector2::Zero();
    float mass = 1.0f;

    for (int otherParticleIndex : particleIndices) {
        if (index == otherParticleIndex) continue;

        Vector2 offset = positions[otherParticleIndex] - positions[index];
        float sqrDst = offset.sqrMagnitude();

        if (sqrDst > sqrSmoothingRadius) continue;

        float dst = std::sqrt(sqrDst);
        Vector2 dir = (dst == 0) ? getRandomDir() : offset / dst;
        float slope = SmoothingKernelDerivative(dst, smoothingRadius);
        float density = densities[otherParticleIndex];
        float sharedPressure = CalculateSharedPressure(density, densities[index]);
        pressureForce += dir * sharedPressure * slope * mass / density;

        if (show_directional_lines == true)
        {
            glBegin(GL_LINES);
            glColor3f(1.0f, 0.0f, 1.0f);
            glVertex2f(positions[index].X, positions[index].Y);
            glVertex2f(positions[otherParticleIndex].X, positions[otherParticleIndex].Y);
            glEnd();
        }
    }
    return pressureForce;
}
static inline Vector2 InteractionForce(const Vector2& inputPos, float radius, float strength, int particleIndex, const std::vector<Vector2>& positions, int mouseButton) {
    Vector2 interactionForce = Vector2::Zero();
    Vector2 offset = inputPos - positions[particleIndex];
    float sqrDst = offset.sqrMagnitude();
    if (sqrDst < radius * radius)
    {
        float dst = std::sqrt(sqrDst);
        Vector2 dirToInput = (dst < std::numeric_limits<float>::epsilon()) ? Vector2::Zero() : offset / dst;
        float centerT = 1 - dst / radius;
        if (mouseButton == GLFW_MOUSE_BUTTON_LEFT) {
            interactionForce += dirToInput * strength;
        }
        else if (mouseButton == GLFW_MOUSE_BUTTON_RIGHT) {
            interactionForce -= dirToInput * strength;
        }
    }
    //std::cout << interactionForce.X << " , " << interactionForce.Y << std::endl;
    return interactionForce;
}
static inline Vector2 CalculateViscosityForce(int particleIndex, const std::vector<int>& neighborIndices, const std::vector<Vector2>& positions, const std::vector<Vector2>& velocities, float smoothingRadius, float viscosityStrength)
{
    Vector2 viscosityForce = Vector2::Zero();
    Vector2 position = positions[particleIndex];

    for (int otherParticleIndex : neighborIndices)
    {
        if (particleIndex == otherParticleIndex) continue;

        Vector2 offset = positions[otherParticleIndex] - position;
        float dst = offset.magnitude();
        float influence = ViscositySmoothingKernel(dst, smoothingRadius);
        viscosityForce += (velocities[otherParticleIndex] - velocities[particleIndex]) * influence;
    }

    return viscosityForce * viscosityStrength;
}

//Updating and Resetting
std::vector<Ball> balls;
static void initializeBalls(float spacingX, float spacingY, int rows, int cols, int radius, float width, float height, std::vector<Ball>& balls, std::vector<Vector2>& positions) {
    balls.clear();
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            float x = (j + 1) * spacingX;
            float y = (i + 1) * spacingY;
            positions.emplace_back(width + x, height + y);
            velocities.emplace_back(Vector2::Zero());
        }
    }
}
static void updateBallPositions(float spacingX, float spacingY, int rows, int cols, int radius, float width, float height, std::vector<Ball>& balls, std::vector<Vector2>& positions) {
    balls.clear();
    initializeBalls(spacingX, spacingY, rows, cols, radius, width, height, balls, positions);
}
static void resetSimulation(std::vector<Vector2>& positions, std::vector<Vector2>& velocities, std::vector<Ball>& balls, int rows, int cols) {
    positions.clear();
    velocities.clear();
    balls.clear();
    positions.reserve(rows * cols);
    velocities.reserve(rows * cols);
    balls.reserve(rows * cols);

    if (reset_to_grid == true) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                float x = (j + 1) * spacingX;
                float y = (i + 1) * spacingY;
                positions.emplace_back(width + x, height + y);
                velocities.emplace_back(Vector2::Zero());
            }
        }
    }
    if (reset_to_random == true) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                positions.emplace_back(distX(mt), distY(mt));
                velocities.emplace_back(Vector2::Zero());
            }
        }
    }

    for (size_t i = 0; i < positions.size(); ++i) {
        balls.emplace_back(positions[i], velocities[i], radius);
    }
}

//Drawing
static void drawDensityAreas(const std::vector<float>& densities, float targetDensity) {
    glBegin(GL_QUADS);
    for (int i = 0; i < densities.size(); ++i) {
        float densityError = densities[i] - targetDensity;
        if (densityError > 0.0001) {
            glColor3f(1.0f, 0.0f, 0.0f); // Red, positive
        }
        else if (densityError < -0.0001) {
            glColor3f(0.0f, 0.0f, 1.0f); // Blue, negative
        }
        else {
            glColor3f(1.0f, 1.0f, 1.0f); // White
        }

        // Calculate the coordinates for the quad based on ball position and radius
        const Vector2& position = balls[i].position;
        float radius = balls[i].radius;
        glVertex2f(position.X - radius, position.Y - radius);
        glVertex2f(position.X + radius, position.Y - radius);
        glVertex2f(position.X + radius, position.Y + radius);
        glVertex2f(position.X - radius, position.Y + radius);
    }
    glEnd();
}
static Vector2 CalculateDensityGradient(const std::vector<Vector2>& positions, float smoothingRadius, int index) {
    const float mass = 1;
    const float epsilon = 1e-6f; // Small value to prevent division by zero

    Vector2 gradient(0.0f, 0.0f);
    Vector2 position = positions[index];

    for (const auto& p : positions) {
        Vector2 offset = p - position;
        float dst = offset.magnitude();

        if (dst < epsilon || dst > smoothingRadius) continue;

        float influence = SmoothingKernel(smoothingRadius, dst);
        Vector2 partialGradient = offset * (mass * influence / (dst + epsilon)); // Avoid division by zero
        gradient += partialGradient;
    }

    gradient.normalize(); // Normalize the final gradient vector
    return gradient;
}
static void DrawDensityGradients(const std::vector<Vector2>& positions, float smoothingRadius) {
    const float arrowLengthFactor = 20.0f; // Adjust this factor for arrow length
    const float arrowHeadSize = 2.0f; // Adjust this size for arrow head

    for (int i = 0; i < positions.size(); ++i) {
        Vector2 gradient = CalculateDensityGradient(positions, smoothingRadius, i);

        // Calculate arrow length based on density gradient magnitude
        float arrowLength = gradient.magnitude() * arrowLengthFactor;

        // Draw the arrow with increased thickness
        glLineWidth(2.0f);
        glBegin(GL_LINES);
        glVertex2f(positions[i].X, positions[i].Y);
        glVertex2f(positions[i].X + gradient.X * arrowLength, positions[i].Y + gradient.Y * arrowLength);
        glEnd();
        glLineWidth(1.0f); // Reset line width

        // Calculate arrowhead vertices
        // Calculate arrowhead vertices
        Vector2 arrowEnd = positions[i] + gradient * arrowLength;
        // Calculate the magnitude of the gradient vector
        float gradientMag = sqrt(gradient.X * gradient.X + gradient.Y * gradient.Y);

        // Calculate the normalized gradient vector components
        float normalizedX = gradient.X / gradientMag;
        float normalizedY = gradient.Y / gradientMag;

        // Create the normalized gradient vector
        Vector2 arrowDir(normalizedX, normalizedY);

        Vector2 arrowHead1 = arrowEnd + Vector2(-arrowDir.Y, arrowDir.X) * arrowHeadSize; // Rotate by 90 degrees
        Vector2 arrowHead2 = arrowEnd + Vector2(arrowDir.Y, -arrowDir.X) * arrowHeadSize; // Rotate by -90 degrees

        // Draw arrow head
        glBegin(GL_TRIANGLES);
        glVertex2f(arrowEnd.X, arrowEnd.Y);
        glVertex2f(arrowHead1.X, arrowHead1.Y);
        glVertex2f(arrowHead2.X, arrowHead2.Y);
        glEnd();

    }
}

//Optimized Particle Lookup
static inline std::pair<float, float> PositionToCellCoord(const std::pair<float, float>& point, float radius) {
    //Convert a position to the coordinate of the cell it is within
    float cellX = static_cast<float>(point.first / radius);
    float cellY = static_cast<float>(point.second / radius);
    return std::make_pair(cellX, cellY);
}
static inline unsigned int HashCell(int cellX, int cellY) {
    // Convert a cell coordinate into a single number.
    unsigned int a = static_cast<unsigned int>(cellX) * 15823u;
    unsigned int b = static_cast<unsigned int>(cellY) * 9737333u;
    return a + b;
}
//static inline unsigned int GetKeyFromHash(unsigned int hash, size_t spatialLookupLength) {
//    // wrap the hash value around the length of the array (so it can be used as an index)
//    return hash % spatialLookupLength;
//}
static inline unsigned int GetKeyFromHash(unsigned int hash, size_t spatialLookupLength) {
    return hash & (spatialLookupLength - 1);
}
static struct Entry {
    int index;
    uint32_t cellKey;
    int cellX;
    int cellY;
    Entry() : index(0), cellKey(0), cellX(0), cellY(0) {}  // Default constructor
    Entry(int index, uint32_t cellKey, int cellX, int cellY) : index(index), cellKey(cellKey), cellX(cellX), cellY(cellY) {}
};

// Forward declaration
static inline void DrawGridCellOutlines(const std::vector<Entry>& spatialLookup, float cellSize);
static inline void HighlightGridCellOutlines(const std::vector<Entry>& spatialLookup, float cellSize, int highlightCellX, int highlightCellY);
static inline void RadixSort(std::vector<Entry>& spatialLookup) {
    const int numBits = 8;
    const int numBuckets = 1 << numBits;
    const int numPasses = sizeof(uint32_t) * 8 / numBits;

    std::vector<Entry> temp(spatialLookup.size());
    std::vector<int> count(numBuckets);

    for (int pass = 0; pass < numPasses; ++pass) {
        int shift = pass * numBits;
        int mask = numBuckets - 1;

        count.assign(numBuckets, 0);

        for (const auto& entry : spatialLookup) {
            int bucket = (entry.cellKey >> shift) & mask;
            ++count[bucket];
        }

        for (int i = 1; i < numBuckets; ++i) {
            count[i] += count[i - 1];
        }

        for (int i = spatialLookup.size() - 1; i >= 0; --i) {
            int bucket = (spatialLookup[i].cellKey >> shift) & mask;
            temp[--count[bucket]] = spatialLookup[i];
        }

        std::swap(spatialLookup, temp);
    }
}
static inline void UpdateSpatialLookup(const std::vector<Vector2>& positions, float radius, std::vector<Entry>& spatialLookup, std::vector<int>& startIndices) {
    spatialLookup.resize(positions.size());

    for (size_t i = 0; i < positions.size(); ++i) {
        int cellX, cellY;
        std::tie(cellX, cellY) = PositionToCellCoord({ positions[i].X, positions[i].Y }, radius);
        uint32_t cellKey = GetKeyFromHash(HashCell(cellX, cellY), spatialLookup.size());
        spatialLookup[i] = Entry(static_cast<int>(i), cellKey, cellX, cellY);
    }

    RadixSort(spatialLookup);

    for (size_t i = 0; i < positions.size(); ++i) {
        uint32_t key = spatialLookup[i].cellKey;
        uint32_t keyPrev = (i == 0) ? UINT32_MAX : spatialLookup[i - 1].cellKey;
        if (key != keyPrev) {
            startIndices[key] = static_cast<int>(i);
        }
    }

    if (show_spatial_grid == true) {
        DrawGridCellOutlines(spatialLookup, smoothingRadius);
    }
}
std::vector<int> particleIndices;
static inline std::array<std::pair<int, int>, 9> cellOffsets = { {
    {-1, -1}, {0, -1}, {1, -1},
    {-1,  0}, {0,  0}, {1,  0},
    {-1,  1}, {0,  1}, {1,  1}
} };
static inline std::vector<int> ForEachPointWithinRadius(const Vector2& samplePoint, float radius, const std::vector<Entry>& spatialLookup, const std::vector<int>& startIndices, const std::vector<Vector2>& positions) {
    particleIndices.clear();
    // Find which cell the sample point is in (this will be the center of our 3x3 block)
    std::pair<int, int> center = PositionToCellCoord({ samplePoint.X, samplePoint.Y }, radius);
    int centerX = center.first;
    int centerY = center.second;

    // Loop over all cells of the 3x3 block around the cell
    for (const auto& offset : cellOffsets) {
        int offsetX = offset.first;
        int offsetY = offset.second;
        // Get key of current cell, then loop over all points that share that key
        uint32_t key = GetKeyFromHash(HashCell(centerX + offsetX, centerY + offsetY), spatialLookup.size());
        int cellStartIndex = startIndices[key];

        //std::cout << "  Cell: (" << centerX + offsetX << ", " << centerY + offsetY << "), key = " << key << ", start index = " << cellStartIndex << std::endl;

        for (size_t i = cellStartIndex; i < spatialLookup.size(); ++i) {
            // Exit loop if we're no longer looking at the correct cell
            if (spatialLookup[i].cellKey != key) break;

            int particleIndex = spatialLookup[i].index;
            float dx = positions[particleIndex].X - samplePoint.X;
            float dy = positions[particleIndex].Y - samplePoint.Y;
            float sqrDst = dx * dx + dy * dy;
            // Test if the point is inside the radius
            if (abs(dx) <= radius && abs(dy) <= radius) {
                float sqrDst = dx * dx + dy * dy;
                if (sqrDst <= sqrSmoothingRadius) {
                    particleIndices.push_back(particleIndex);
                }
            }
        }
    }
    //std::cout << "  Found " << particleIndices.size() << " particles within radius" << std::endl;
    return particleIndices;
}

//Drawing
static void DrawGridCellOutlines(const std::vector<Entry>& spatialLookup, float cellSize) {
    std::set<std::pair<int, int>> drawnCells;

    for (const auto& entry : spatialLookup) {
        if (drawnCells.find({ entry.cellX, entry.cellY }) == drawnCells.end()) {
            // This cell has not been drawn yet, so draw it now.
            drawnCells.insert({ entry.cellX, entry.cellY });

            // Calculate the bottom left corner of the cell
            float x = entry.cellX * cellSize;
            float y = entry.cellY * cellSize;


            glColor3f(0.5f, 0.5f, 0.5f);
            glBegin(GL_LINE_LOOP);
            glVertex2f(x, y);
            glVertex2f(x + cellSize, y);
            glVertex2f(x + cellSize, y + cellSize);
            glVertex2f(x, y + cellSize);
            glEnd();
        }
    }
}
static void HighlightGridCellOutlines(const std::vector<Entry>& spatialLookup, float cellSize, int highlightCellX, int highlightCellY) {
    std::set<std::pair<int, int>> drawnCells;

    for (const auto& entry : spatialLookup) {
        if (entry.cellX == highlightCellX && entry.cellY == highlightCellY) {
            // Calculate the bottom left corner of the cell
            float x = entry.cellX * cellSize;
            float y = entry.cellY * cellSize;

            glColor3f(1.0f, 0.5f, 0.0f);
            glBegin(GL_LINE_LOOP);
            glVertex2f(x, y);
            glVertex2f(x + cellSize, y);
            glVertex2f(x + cellSize, y + cellSize);
            glVertex2f(x, y + cellSize);
            glEnd();

            break;
        }
    }
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

    //glfwSwapInterval(1);

    std::vector<Vector2> positions;
    std::vector<Vector2> predictedPositions;
    positions.reserve(rows * cols);
    std::vector<Vector2> velocities(rows * cols, Vector2(0.0, 0.0));
    std::vector<Vector2> mousePosition;
    /*for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            positions.emplace_back(distX(mt), distY(mt));
        }
    }*/
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            float x = (j + 1) * spacingX;
            float y = (i + 1) * spacingY;
            positions.emplace_back(width + x, height + y);
        }
    }
    balls.reserve(rows * cols);
    for (size_t i = 0; i < positions.size(); ++i) {
        balls.emplace_back(positions[i], velocities[i], radius);
    }
    
    std::vector<float> densities;
    std::vector<Entry> spatialLookup(positions.size());
    std::vector<int> startIndices(positions.size(), std::numeric_limits<int>::max());
    predictedPositions = positions;
    std::vector<bool> isParticleLost(balls.size(), false);
    std::vector<std::vector<int>> particleIndices(positions.size());
    densities.resize(positions.size());

    std::cout << "Entering window loop" << std::endl;
    while (!glfwWindowShouldClose(window))
    {
        glClear(GL_COLOR_BUFFER_BIT);

        glfwGetCursorPos(window, &mouseX, &mouseY);
        Vector2 mousePosition(static_cast<float>(mouseX), static_cast<float>(-mouseY + SCREEN_HEIGHT));
        if (show_mouse == true) {
            drawOutline(mousePosition, mouseRadius);
        }

        for (int i = 0; i < balls.size(); i++) {
            velocities[i] += Vector2::down() * gravity;
            predictedPositions[i] = positions[i] + velocities[i];
        }
        UpdateSpatialLookup(predictedPositions, smoothingRadius, spatialLookup, startIndices);

        float density = CalculateDensity(positions, smoothingRadius);
        for (int i = 0; i < positions.size(); ++i) {
            particleIndices[i] = ForEachPointWithinRadius(predictedPositions[i], smoothingRadius, spatialLookup, startIndices, predictedPositions);
        }

        for (int i = 0; i < positions.size(); ++i) {
            densities[i] = CalculateDensity(i, particleIndices[i], predictedPositions, smoothingRadius);
        }

        std::for_each(std::execution::par, positions.begin(), positions.end(), [&](const auto& position) {
            int i = &position - &positions[0];
            Vector2 pressureForce = CalculatePressureForce(i, particleIndices[i], predictedPositions, densities, smoothingRadius);
            Vector2 pressureAcceleration = pressureForce / densities[i];
            velocities[i] += pressureAcceleration;
            });

        std::for_each(std::execution::par, positions.begin(), positions.end(), [&](const auto& position) {
            int i = &position - &positions[0];
            Vector2 viscosityForce = CalculateViscosityForce(i, particleIndices[i], predictedPositions, velocities, smoothingRadius, viscosityStrength);
            velocities[i] += viscosityForce;
            });

        int mouseButton = -1;
        if (enableInteraction == true)
        {
            if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
                mouseButton = GLFW_MOUSE_BUTTON_LEFT;
            }
            else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
                mouseButton = GLFW_MOUSE_BUTTON_RIGHT;
            }
        }

        std::for_each(std::execution::par, positions.begin(), positions.end(), [&](const auto& position) {
            int i = &position - &positions[0];
            float density = std::max(densities[i], 1.0f);
            Vector2 interactionForce = InteractionForce(mousePosition, mouseRadius, mouseStrength, i, predictedPositions, mouseButton);
            velocities[i] += interactionForce / density;
            });

        std::for_each(std::execution::par, positions.begin(), positions.end(), [&](const auto& position) {
            int i = &position - &positions[0];
            positions[i] += velocities[i];
            balls[i].position = positions[i];
            balls[i].velocity = velocities[i];
            });
    
        //for (int i = 0; i < balls.size(); i++) {
        //    if (!isParticleLost[i] && (std::isnan(velocities[i].X) || std::isnan(velocities[i].Y))) {
        //        lostParticles++;
        //        isParticleLost[i] = true; // This particle is now counted as lost
        //    }
        //}

        // Define obstacle properties
        
        // Resolve collisions with the bounding box and obstacle
        resolveCollisions(positions, velocities, radius, dampingFactor, boundsSizeX, boundsSizeY, obstacleCenter, obstacleWidth, obstacleHeight);

        // ...

        for (const Ball& ball : balls) {
            ball.draw();
        }
        //DrawDensityGradients(positions, smoothingRadius);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        {
            static int ballRadius = radius;
            ImGui::Begin("Config");
            ImGui::SliderFloat("Gravity", &gravity, 0.0f, 0.2f);
            ImGui::SliderInt("Ball Radius", &ballRadius, 1, 70);
            ImGui::Text("Density: %.8f", density);
            ImGui::SliderFloat("Smoothing Radius", &smoothingRadius, 60.0f, 500.0f);
            //ImGui::SliderInt("Spacing X", &spacingX, 0.0, 150.0);
            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
            ImGui::SliderFloat("Pressure Multiplier", &pressureMultiplier, 50.0f, 1000.0f);
            //ImGui::SliderFloat("Multiplicative Factor", &multiplicativeFactor, 1.0f, 100.0f);
            ImGui::SliderFloat("Target Density", &targetDensity, 0.000001f, 0.01f, "%.8f");
            ImGui::SliderFloat("Mouse Radius", &mouseRadius, 10.0f, 200.0f);
            ImGui::SliderFloat("Mouse Strength", &mouseStrength, 0.1f, 2.0f);
            ImGui::SliderFloat("Viscosity Strength", &viscosityStrength, 0.1f, 20.0f);
            ImGui::SliderFloat2("Obstacle Center", &obstacleCenter.X, 0.0f, boundsSizeX);
            ImGui::SliderFloat("Obstacle Width", &obstacleWidth, 10.0f, boundsSizeX);
            ImGui::SliderFloat("Obstacle Height", &obstacleHeight, 10.0f, boundsSizeY);
            ImGui::Checkbox("Show Smoothing Radius", &show_smoothing_radius);
            ImGui::SameLine();
            ImGui::Checkbox("Show Directional Lines", &show_directional_lines);
            ImGui::Checkbox("Show Density Areas", &show_density_areas);
            ImGui::SameLine();
            ImGui::Checkbox("Show Spatial Grid", &show_spatial_grid);
            ImGui::Checkbox("Enable Interaction", &enableInteraction);

            //ImGui::Text("Lost Particles: %d", lostParticles);
            if (ImGui::Button("Reset to Grid")) {
                reset_to_grid = true;
                reset_to_random = false;
                resetSimulation(positions, velocities, balls, rows, cols);
            }
            ImGui::SameLine();
            if (ImGui::Button("Reset to Random")) {
                reset_to_grid = false;
                reset_to_random = true;
                resetSimulation(positions, velocities, balls, rows, cols);
            }
            ImGui::End();

            if (radius != ballRadius) {
                radius = ballRadius;
                for (auto& ball : balls) ball.radius = radius;
            }
            if (show_smoothing_radius) drawBounds(smoothingRadius);
            if (show_directional_lines) show_directional_lines = true;
            if (enableInteraction) enableInteraction = true;
            if (show_density_areas) { drawDensityAreas(densities, targetDensity); show_density_areas = true; }
        }
        if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {
            resetSimulation(positions, velocities, balls, rows, cols);
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

//{
//    Vector2 mousePosition(static_cast<float>(mouseX), static_cast<float>(-mouseY + SCREEN_HEIGHT));
//    if (show_mouse == true) {
//        drawOutline(mousePosition, mouseRadius);
//    }
//
//    for (int i = 0; i < balls.size(); i++) {
//        velocities[i] += Vector2::down() * gravity;
//        predictedPositions[i] = positions[i] + velocities[i];
//    }
//
//    UpdateSpatialLookup(predictedPositions, smoothingRadius, spatialLookup, startIndices);
//
//    float density = CalculateDensity(positions, smoothingRadius);
//    for (int i = 0; i < positions.size(); ++i) {
//        particleIndices[i] = ForEachPointWithinRadius(predictedPositions[i], smoothingRadius, spatialLookup, startIndices, predictedPositions);
//    }
//
//    for (int i = 0; i < positions.size(); ++i) {
//        densities[i] = CalculateDensity(i, particleIndices[i], predictedPositions, smoothingRadius);
//    }
//
//    /*for (int i = 0; i < positions.size(); ++i) {
//
//        Vector2 pressureForce = CalculatePressureForce(i, particleIndices[i], predictedPositions, densities, smoothingRadius);
//        Vector2 pressureAcceleration = pressureForce / densities[i];
//        velocities[i] += pressureAcceleration;
//    }*/
//
//    std::for_each(std::execution::par, positions.begin(), positions.end(), [&](const auto& position) {
//        int i = &position - &positions[0];
//        Vector2 pressureForce = CalculatePressureForce(i, particleIndices[i], predictedPositions, densities, smoothingRadius);
//        Vector2 pressureAcceleration = pressureForce / densities[i];
//        velocities[i] += pressureAcceleration;
//        });
//
//    for (int i = 0; i < positions.size(); ++i) {
//
//        Vector2 viscosityForce = CalculateViscosityForce(i, particleIndices[i], predictedPositions, velocities, smoothingRadius, viscosityStrength);
//        velocities[i] += viscosityForce;
//    }
//
//    int mouseButton = -1;
//    if (enableInteraction == true)
//    {
//        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
//            mouseButton = GLFW_MOUSE_BUTTON_LEFT;
//        }
//        else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
//            mouseButton = GLFW_MOUSE_BUTTON_RIGHT;
//        }
//    }
//
//    for (int i = 0; i < balls.size(); i++) {
//        float density = std::max(densities[i], 1.0f); // Clamp density to a minimum value of 0.01
//        Vector2 interactionForce = InteractionForce(mousePosition, mouseRadius, mouseStrength, i, predictedPositions, mouseButton);
//        velocities[i] += interactionForce / density;
//    }
//
//    for (int i = 0; i < balls.size(); i++) {
//        positions[i] += velocities[i];
//        balls[i].position = positions[i];
//        balls[i].velocity = velocities[i]; // Update the velocity of each Ball object
//    }
//    //for (int i = 0; i < balls.size(); i++) {
//    //    if (!isParticleLost[i] && (std::isnan(velocities[i].X) || std::isnan(velocities[i].Y))) {
//    //        lostParticles++;
//    //        isParticleLost[i] = true; // This particle is now counted as lost
//    //    }
//    //}
//
//    resolveCollisions(positions, velocities, radius, dampingFactor, boundsSizeX, boundsSizeY);
//
//    for (const Ball& ball : balls) {
//        ball.draw();
//    }
//}