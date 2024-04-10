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
#ifndef M_PI
#define M_PI 3.14
#endif
GLfloat twicePi = 2.0f * M_PI;
constexpr auto SCREEN_WIDTH = 1280;
constexpr auto SCREEN_HEIGHT = 720;

int rows = 20;
int cols = rows;
int radius = 4;
int spacingX = 11;
int spacingY = spacingX;
float smoothingRadius = 66.0f;
float gravity = 0.000f;
float dampingFactor = 0.8f;
float targetDensity = 0.0003f;
float pressureMultiplier = 2.80f;
float stiffnessConstant = 1.0f;
float boundsSizeX = SCREEN_WIDTH;
float boundsSizeY = SCREEN_HEIGHT;
float width = SCREEN_WIDTH / 2.0 - (((cols + 1.0) * spacingX) / 2.0);
float height = SCREEN_HEIGHT / 2.0 - (((cols + 1.0) * spacingX) / 2.0);
bool show_smoothing_radius = false;
bool show_directional_lines = false;
bool show_density_areas = false;
bool show_spatial_grid = false;
float mouseRadius = smoothingRadius;
GLint numberOfSides = 32;
std::mt19937 mt{ std::random_device{}() };
std::uniform_real_distribution<float> distX(0.0f, SCREEN_WIDTH);
std::uniform_real_distribution<float> distY(0.0f, SCREEN_HEIGHT);
float sqrSmoothingRadius = smoothingRadius * smoothingRadius;

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

static float SmoothingKernel(float radius, float dst)
{
    if (dst < radius)
    {
        float volume = (M_PI * pow(radius, 4)) / 6;
        float result = ((radius - dst) * (radius - dst) / volume);
        
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
std::vector<Vector2> velocities;

class Ball {
public:
    Vector2 position;
    Vector2 velocity;
    float radius;

    Ball(const Vector2& position, const Vector2& velocity, float radius)
        : position(position), velocity(velocity), radius(radius) {}

    void draw() const {
        drawCircle(position.X, position.Y, 0, radius, 0.0, 0.8, 1.0);
    }
    void draw(GLfloat r, GLfloat g, GLfloat b) const {
        drawCircle(position.X, position.Y, 0, radius, r, g, b);
    }
    void drawOutlineCircle() const {
        drawCircle(position.X, position.Y, 0, radius, 0.0, 1.0, 0.0);
    }

public:
    void drawCircle(GLfloat x, GLfloat y, GLfloat z, GLfloat radius, GLfloat r, GLfloat g, GLfloat b) const {
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

void resolveCollisions(std::vector<Vector2>& positions, std::vector<Vector2>& velocities, float radius, float dampingFactor, float boundsSizeX, float boundsSizeY) {
    for (int i = 0; i < positions.size(); i++) {
        float minX = radius;
        float maxX = boundsSizeX - radius;
        float minY = radius;
        float maxY = boundsSizeY - radius;

        if (positions[i].X < minX) {
            positions[i].X = minX;
            velocities[i].X *= -dampingFactor;
        }
        else if (positions[i].X > maxX) {
            positions[i].X = maxX;
            velocities[i].X *= -dampingFactor;
        }

        if (positions[i].Y < minY) {
            positions[i].Y = minY;
            velocities[i].Y *= -dampingFactor;
        }
        else if (positions[i].Y > maxY) {
            positions[i].Y = maxY;
            velocities[i].Y *= -dampingFactor;
        }
    }
}

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
void drawOutline(Vector2 position, float smoothingRadius) {
    glColor3f(0.0f, 1.0f, 0.0f);
    //drawRadius(position.X, position.Y, 0, radius, 30);
    drawRadius(position.X, position.Y, 0, smoothingRadius, 128);
}

float CalculateDensity(const std::vector<Vector2>& positions, float smoothingRadius, int index) {
    float density = 0;
    const float mass = 1;
    Vector2 samplePoint(SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2);
    Vector2 position = positions[index];
    for (const auto& p : positions) {
        Vector2 offset = p - position;
        float sqrDst = offset.sqrMagnitude();

        if (sqrDst > sqrSmoothingRadius) continue;

        float influence = SmoothingKernel(smoothingRadius, std::sqrt(sqrDst));
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
    return (density - targetDensity) * pressureMultiplier;
}
static float CalculateDensityError(float density)
{
    float densityError = density - targetDensity;
    return densityError;

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
    float mass = 1.0f;

    for (int otherParticleIndex = 0; otherParticleIndex < positions.size(); otherParticleIndex++)
    {
        if (particleIndex == otherParticleIndex) continue;

        Vector2 offset = positions[otherParticleIndex] - positions[particleIndex];
        float sqrDst = offset.sqrMagnitude();

        if (sqrDst > sqrSmoothingRadius) continue;

        float dst = std::sqrt(sqrDst);
        Vector2 dir = (dst == 0) ? getRandomDir() : offset / dst;
        float slope = SmoothingKernelDerivative(dst, smoothingRadius);
        float density = densities[otherParticleIndex];
        float sharedPressure = CalculateSharedPressure(density, densities[particleIndex]);
        pressureForce += dir * sharedPressure * slope * mass / density;

        if (show_directional_lines == true)
        {
            glBegin(GL_LINES);
            glColor3f(1.0f, 0.0f, 1.0f);
            glVertex2f(positions[particleIndex].X, positions[particleIndex].Y);
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

//void initializeBalls(float spacingX, float spacingY, int rows, int cols, int radius, float width, float height, std::vector<Ball>& balls) {
//    balls.clear();
//    for (int i = 0; i < rows; ++i) {
//        for (int j = 0; j < cols; ++j) {
//            float x = (j + 1) * spacingX;
//            float y = (i + 1) * spacingY;
//            balls.emplace_back(width + x, height + y, 0.0, 0.0, radius);
//        }
//    }
//}

//void updateBallPositions(float spacingX, float spacingY, int rows, int cols, int radius, float width, float height, std::vector<Ball>& balls) {
//    balls.clear();
//    initializeBalls(spacingX, spacingY, rows, cols, radius, width, height, balls);
//}
void resetSimulation(std::vector<Vector2>& positions, std::vector<Vector2>& velocities, std::vector<Ball>& balls, int rows, int cols) {
    positions.clear();
    velocities.clear();
    balls.clear();
    positions.reserve(rows * cols);
    velocities.reserve(rows * cols);
    balls.reserve(rows * cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            positions.emplace_back(distX(mt), distY(mt));
            velocities.emplace_back(Vector2::Zero());
        }
    }
    /*for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            float x = (j + 1) * spacingX;
            float y = (i + 1) * spacingY;
            positions.emplace_back(width + x, height + y);
            velocities.emplace_back(Vector2::Zero());
        }
    }*/

    for (size_t i = 0; i < positions.size(); ++i) {
        balls.emplace_back(positions[i], velocities[i], radius);
    }
}

std::vector<Ball> balls;
void drawDensityAreas(const std::vector<float>& densities, float targetDensity) {
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

Vector2 CalculateDensityGradient(const std::vector<Vector2>& positions, float smoothingRadius, int index) {
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

void DrawDensityGradients(const std::vector<Vector2>& positions, float smoothingRadius) {
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


//Spatial Grid
//Convert a position to the coordinate of the cell it is within
std::pair<float, float> PositionToCellCoord(const std::pair<float, float>& point, float radius) {
    float cellX = static_cast<float>(point.first / radius);
    float cellY = static_cast<float>(point.second / radius);
    return std::make_pair(cellX, cellY);
}


// Convert a cell coordinate into a single number.
unsigned int HashCell(int cellX, int cellY) {
    unsigned int a = static_cast<unsigned int>(cellX) * 15823u;
    unsigned int b = static_cast<unsigned int>(cellY) * 9737333u;
    return a + b;
}

// wrap the hash value around the length of the array (so it can be used as an index)
unsigned int GetKeyFromHash(unsigned int hash, size_t spatialLookupLength) {
    return hash % spatialLookupLength;
}

struct Entry {
    int index;
    uint32_t cellKey;
    int cellX;
    int cellY;

    Entry() : index(0), cellKey(0), cellX(0), cellY(0) {}  // Default constructor

    Entry(int index, uint32_t cellKey, int cellX, int cellY) : index(index), cellKey(cellKey), cellX(cellX), cellY(cellY) {}
};

// Forward declaration
void DrawGridCellOutlines(const std::vector<Entry>& spatialLookup, float cellSize);
void HighlightGridCellOutlines(const std::vector<Entry>& spatialLookup, float cellSize, int highlightCellX, int highlightCellY);

void UpdateSpatialLookup(const std::vector<Vector2>& positions, float radius, std::vector<Entry>& spatialLookup, std::vector<int>& startIndices) {
    // Create (unordered) spatial lookup
    spatialLookup.resize(positions.size());

    //std::cout << "Updating spatial lookup with " << positions.size() << " particles" << std::endl;

    for (size_t i = 0; i < positions.size(); ++i) {
        int cellX, cellY;
        std::tie(cellX, cellY) = PositionToCellCoord({ positions[i].X, positions[i].Y }, radius);
        uint32_t cellKey = GetKeyFromHash(HashCell(cellX, cellY), spatialLookup.size());
        spatialLookup[i] = Entry(static_cast<int>(i), cellKey, cellX, cellY);
        //std::cout << "  Particle " << i << ": cell = (" << cellX << ", " << cellY << "), key = " << cellKey << std::endl;
    }

    std::sort(spatialLookup.begin(), spatialLookup.end(), [](const Entry& a, const Entry& b) {
        return a.cellKey < b.cellKey;
        });

    // Calculate start indices of each unique cell key in the spatial lookup
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

std::vector<int> ForEachPointWithinRadius(const Vector2& samplePoint, float radius, const std::vector<Entry>& spatialLookup, const std::vector<int>& startIndices, const std::vector<Vector2>& positions) {
    //std::cout << "ForEachPointWithinRadius: samplePoint = (" << samplePoint.X << ", " << samplePoint.Y << "), radius = " << radius << std::endl;

    // Find which cell the sample point is in (this will be the center of our 3x3 block)
    std::pair<int, int> center = PositionToCellCoord({ samplePoint.X, samplePoint.Y }, radius);
    int centerX = center.first;
    int centerY = center.second;
    float sqrRadius = radius * radius;
    std::vector<int> particleIndices;

    // cellOffsets array
    std::array<std::pair<int, int>, 9> cellOffsets = { {
        {-1, -1}, {0, -1}, {1, -1},
        {-1,  0}, {0,  0}, {1,  0},
        {-1,  1}, {0,  1}, {1,  1}
    } };

    // Loop over all cells of the 3x3 block around the cell
    for (const auto& offset : cellOffsets) {
        int offsetX = offset.first;
        int offsetY = offset.second;
        // Get key of current cell, then loop over all points that share that key
        uint32_t key = GetKeyFromHash(HashCell(centerX + offsetX, centerY + offsetY), spatialLookup.size());
        int cellStartIndex = startIndices[key];

        //std::cout << "  Cell: (" << centerX + offsetX << ", " << centerY + offsetY << "), key = " << key << ", start index = " << cellStartIndex << std::endl;
        
        if (show_spatial_grid == true) {
            HighlightGridCellOutlines(spatialLookup, smoothingRadius, centerX + offsetX, centerY + offsetY);
        }
        for (size_t i = cellStartIndex; i < spatialLookup.size(); ++i) {
            // Exit loop if we're no longer looking at the correct cell
            if (spatialLookup[i].cellKey != key) break;

            int particleIndex = spatialLookup[i].index;
            float dx = positions[particleIndex].X - samplePoint.X;
            float dy = positions[particleIndex].Y - samplePoint.Y;
            float sqrDst = dx * dx + dy * dy;
            // Test if the point is inside the radius
            if (sqrDst <= sqrRadius) {
                particleIndices.push_back(particleIndex);
                //std::cout << "    Particle " << particleIndex << " is within radius: (" << positions[particleIndex].X << ", " << positions[particleIndex].Y << ")" << std::endl;
            }
        }
    }

    //std::cout << "  Found " << particleIndices.size() << " particles within radius" << std::endl;
    return particleIndices;
}


void DrawGridCellOutlines(const std::vector<Entry>& spatialLookup, float cellSize) {
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

void HighlightGridCellOutlines(const std::vector<Entry>& spatialLookup, float cellSize, int highlightCellX, int highlightCellY) {
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
    positions.reserve(rows * cols); 
    std::vector<Vector2> velocities(rows * cols, Vector2(0.0, 0.0));
    std::vector<Vector2> mousePosition;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            positions.emplace_back(distX(mt), distY(mt));
        }
    }
    /*for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            float x = (j + 1) * spacingX;
            float y = (i + 1) * spacingY;
            positions.emplace_back(width + x, height + y);
        }
    }*/
    balls.reserve(rows * cols);
    for (size_t i = 0; i < positions.size(); ++i) {
        balls.emplace_back(positions[i], velocities[i], radius);
    }

    /*for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            float x = (j + 1) * spacingX;
            float y = (i + 1) * spacingY;
            balls.emplace_back(width + x, height + y, 0.0, 0.0, radius);
        }
    }*/
    std::cout << "Entering window loop" << std::endl;
    std::vector<float> densities;
    positions.reserve(balls.size());
    densities.reserve(balls.size());
    float lastFrame = 0.0f;
    double mouseX, mouseY;
    std::vector<Entry> spatialLookup(positions.size()); // Make sure it's resized with the correct size
    std::vector<int> startIndices(positions.size(), std::numeric_limits<int>::max()); // Similar for startIndices

    while (!glfwWindowShouldClose(window))
    {
        glClear(GL_COLOR_BUFFER_BIT);

        float currentFrame = glfwGetTime();
        
        float deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
        
        glfwGetCursorPos(window, &mouseX, &mouseY);
        Vector2 mousePosition(static_cast<float>(mouseX), static_cast<float>(-mouseY + 720));

        if (show_spatial_grid == true) {
            drawOutline(mousePosition, smoothingRadius);
        }

        UpdateSpatialLookup(positions, smoothingRadius, spatialLookup, startIndices);

        std::vector<int> particleIndices = ForEachPointWithinRadius(mousePosition, smoothingRadius, spatialLookup, startIndices, positions);

        for (int i = 0; i < balls.size(); i++) {
            if (std::find(particleIndices.begin(), particleIndices.end(), i) != particleIndices.end()) {

                if (show_spatial_grid == true) {
                    balls[i].draw(1.0, 0.0, 0.0);  // Red color
                }
                else {
                    balls[i].draw(0.0, 0.8, 1.0);  // Default color
                }
            }
            else {
                balls[i].draw(0.0, 0.8, 1.0);  // Default color
            }
        }

        for (int i = 0; i < balls.size(); i++) {
            velocities[i] += Vector2::down() * gravity;
            //predictedositions[i] = positions[i] * velocities[i];
        }

        float density = CalculateDensity(positions, smoothingRadius);
        densities.resize(positions.size());
        PreCalculateDensities(densities, positions, smoothingRadius);

        for (int i = 0; i < positions.size(); ++i) {
            Vector2 pressureForce = CalculatePressureForce(i, positions, densities, smoothingRadius);
            Vector2 pressureAcceleration = pressureForce / densities[i];
            velocities[i] += pressureAcceleration * stiffnessConstant;
        }

        for (int i = 0; i < balls.size(); i++) {
            positions[i] += velocities[i];
            balls[i].position = positions[i];
        }
        resolveCollisions(positions, velocities, radius, dampingFactor, boundsSizeX, boundsSizeY);

        /*for (auto& ball : balls) {
            ball.draw();
        }*/
        
        //DrawDensityGradients(positions, smoothingRadius);

        //std::cout << deltaTime << std::endl;
        //width = SCREEN_WIDTH / 2 - (((cols + 1) * spacingX) / 2);
        //height = SCREEN_HEIGHT / 2 - (((cols + 1) * spacingY) / 2);
        //updateBallPositions(spacingX, spacingY, rows, cols, radius, width, height, balls);
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        {
            static int ballRadius = radius;
            ImGui::Begin("Config");
            ImGui::SliderFloat("Gravity", &gravity, -0.01f, 0.01f);
            ImGui::SliderInt("Ball Radius", &ballRadius, 1, 70);
            ImGui::Text("Density: %.8f", density);
            ImGui::SliderFloat("Smoothing Radius", &smoothingRadius, 0.05f, 500.0f);
            //ImGui::SliderInt("Spacing X", &spacingX, 0.0, 150.0);
            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
            ImGui::SliderFloat("Pressure Multiplier", &pressureMultiplier, 0.01f, 100.0f);
            //ImGui::SliderFloat("Multiplicative Factor", &multiplicativeFactor, 1.0f, 100.0f);
            ImGui::SliderFloat("Target Density", &targetDensity, 0.000001f, 0.001f, "%.8f");
            ImGui::SliderFloat("Stiffness Constant", &stiffnessConstant, 0.1f, 1.0f);
            ImGui::Checkbox("Show Smoothing Radius", &show_smoothing_radius);
            ImGui::SameLine();
            ImGui::Checkbox("Show Directional Lines", &show_directional_lines);
            ImGui::Checkbox("Show Density Areas", &show_density_areas);
            ImGui::SameLine();
            ImGui::Checkbox("Show Spatial Grid", &show_spatial_grid);
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
            if (show_smoothing_radius) {
                for (const auto& pos : positions) {
                    drawBounds(pos, smoothingRadius);
                }
            }
            if (show_directional_lines) {
                show_directional_lines = true;

            }
            if (show_density_areas) {
                drawDensityAreas(densities, targetDensity);
                show_density_areas = true;
            }
            if (show_density_areas) {
                show_density_areas = true;
            }
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

static void drawLine(float x1, float y1, float x2, float y2) {
    glBegin(GL_LINES);
    glVertex2f(x1, y1);
    glVertex2f(x2, y2);
    glEnd();
}
