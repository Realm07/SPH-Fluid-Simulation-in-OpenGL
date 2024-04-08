#include <tuple>
#include <algorithm>
#include <limits>
#include <utility>
#include <cstdint>

//convert a position to the coordinate of the cell it is within
public (int x, int y) PositionToCellCoord(Vector2 point, float radius)
{
	int cellX = (int)(point.X / radius);
	int cellY = (int)(point.y / radius);
	return (cellX, cellY);
}


//convert a position to the coordinate of the cell it is within
std::tuple<int, int> PositionToCellCoord(Vector2 point, float radius) {
	int cellX = static_cast<int>(point.X / radius);
	int cellY = static_cast<int>(point.Y / radius);
	return std::make_tuple(cellX, cellY);
}






//Convert a cell coordinate into a single number.
public uint HashCell(int cellX, int cellY)
{
	uint a = (uint)cellX * 15823;
	uint b = (uint)cellY * 9737333;
	return a + b;
}

unsigned int HashCell(int cellX, int cellY) {
	unsigned int a = static_cast<unsigned int>(cellX) * 15823u;
	unsigned int b = static_cast<unsigned int>(cellY) * 9737333u;
	return a + b;
}





//wrap the hash value around the length of the array (so it can be used as an index)
public uint GetKeyFromHash(uint hash)
{
	return hash % (uint)spatialLookup.Length;
}

unsigned int GetKeyFromHash(unsigned int hash) {
	unsigned int spatialLookupLength = sizeof(spatialLookup) / sizeof(spatialLookup[0]);
	return hash % spatialLookupLength;
}




public void UpdateSpatialLookup(Vector2[] points, float radius)
{
	this.points = points;
	this.radius = radius;

	//Create (unordered) spatial lookup
	for (0, points.Length, i = >
	{
		(int cellX, int cellY) = PositionToCellCoord(points[i], radius);
		uint cellKey = GetKeyFromHash(HashCell(cellX, cellY));
		spatialLookup[i] = new Entry(i, cellKey);
		startIndices[i] = int.MaxValue; //Reset start index
	});

	Array.Sort(spatialLookup);

	//Calculate start indices of each unique cell key in the spatial lookup
	for (0, points.Length, i = >
	{
		uint key = spatialLookup[i].cellKey;
		uint keyPrev = i == 0 ? uint.MaxValue : spatialLookup[i - 1].cellKey;
		{
			startIndices[key] = i;
		}
	});
}

#include <vector>
#include <algorithm>
#include <limits>
#include <utility>
#include <cstdint>

struct Entry {
	int index;
	uint32_t cellKey;

	Entry(int index, uint32_t cellKey) : index(index), cellKey(cellKey) {}
};

std::pair<int, int> PositionToCellCoord(const std::pair<float, float>& point, float radius);
uint32_t GetKeyFromHash(int hash);
int HashCell(int cellX, int cellY);

void UpdateSpatialLookup(const std::vector<std::pair<float, float>>& points, float radius, std::vector<Entry>& spatialLookup, std::vector<int>& startIndices) {
	// Create (unordered) spatial lookup
	spatialLookup.resize(points.size());
	startIndices.resize(points.size(), std::numeric_limits<int>::max());
	for (size_t i = 0; i < points.size(); ++i) {
		auto [cellX, cellY] = PositionToCellCoord(points[i], radius);
		uint32_t cellKey = GetKeyFromHash(HashCell(cellX, cellY));
		spatialLookup[i] = Entry(i, cellKey);
	}

	std::sort(spatialLookup.begin(), spatialLookup.end(), [](const Entry& a, const Entry& b) {
		return a.cellKey < b.cellKey;
		});

	// Calculate start indices of each unique cell key in the spatial lookup
	for (size_t i = 0; i < points.size(); ++i) {
		uint32_t key = spatialLookup[i].cellKey;
		uint32_t keyPrev = i == 0 ? std::numeric_limits<uint32_t>::max() : spatialLookup[i - 1].cellKey;

		if (key != keyPrev) {
			startIndices[key] = i;
		}
	}
}











public void ForEachPointWithinRadius(Vector2 samplePoint)
{
	//Find which cell the sample point is in (this will be the centre of our 3x3 block)
	float centerX, int centerY = PositionToCellCoord(samplePoint, radius);
	float sqrRadius = radius * radius;

	//Loop over all cells of the 3x3 block around the cell
	foreach((int offsetX, int offsetY) in cellOffsets)
	{
		//Get key of current cell, then loop over all points that share that key
		uint key = GetKeyFromHash(HashCell(centerX + offsetX, centerY + offsetY));
		int cellStartIndex = startIndices[key];

		for (int i = cellStartIndex; i < spatialLookup.Length; i++)
		{
			//exit loop if we're no longer looking at the correct cell
			if (spatialLookup[i].cellKey != key) break;

			int particleIndex = spatialLookup[i].particleIndex;
			float sqrDst = (points[particleIndex] - samplePoint).sqrMagnitude;

			//Test if the point is inside the radius
			if (sqrDst <= sqrRadius)
			{
				//Do something with the particleIndex!
				//either by writing code that uses it directly, or more likely by
				//having this function take in a callback, or return an IEnumerable, etc.
			}
		}
	}
}

#include <vector>
#include <utility>
#include <cstdint>
#include <array>

std::pair<int, int> PositionToCellCoord(const std::pair<float, float>& point, float radius);
uint32_t GetKeyFromHash(int hash);
int HashCell(int cellX, int cellY);

void ForEachPointWithinRadius(const std::pair<float, float>& samplePoint, float radius, const std::vector<Entry>& spatialLookup, const std::vector<int>& startIndices, const std::vector<std::pair<float, float>>& points) {
	// Find which cell the sample point is in (this will be the center of our 3x3 block)
	auto [centerX, centerY] = PositionToCellCoord(samplePoint, radius);
	float sqrRadius = radius * radius;

	// cellOffsets array
	std::array<std::pair<int, int>, 9> cellOffsets = { {
		{-1, -1}, {0, -1}, {1, -1},
		{-1,  0}, {0,  0}, {1,  0},
		{-1,  1}, {0,  1}, {1,  1}
	} };

	// Loop over all cells of the 3x3 block around the cell
	for (const auto& [offsetX, offsetY] : cellOffsets) {
		// Get key of current cell, then loop over all points that share that key
		uint32_t key = GetKeyFromHash(HashCell(centerX + offsetX, centerY + offsetY));
		int cellStartIndex = startIndices[key];

		for (size_t i = cellStartIndex; i < spatialLookup.size(); ++i) {
			// Exit loop if we're no longer looking at the correct cell
			if (spatialLookup[i].cellKey != key) break;

			int particleIndex = spatialLookup[i].index;
			float dx = points[particleIndex].first - samplePoint.first;
			float dy = points[particleIndex].second - samplePoint.second;
			float sqrDst = dx * dx + dy * dy;

			// Test if the point is inside the radius
			if (sqrDst <= sqrRadius) {
				// Do something with the particleIndex!
				// Either by writing code that uses it directly or by
				// having this function take in a callback, or return an IEnumerable, etc.
			}
		}
	}
}

This function takes a std::pair<float, float> as the samplePoint, a float radius, the spatialLookup, startIndices, and points vectors as input.Make sure to implement the helper functions PositionToCellCoord, GetKeyFromHash, and HashCell, as they were not provided in the original code snippet.

In the "// Do something with the particleIndex!" section, you can add the code that you want to execute for each particle index within the radius.













To implement spatial grid optimization in your CalculateDensity function, you can use the UpdateSpatialLookup function provided earlier.First, you need to update the spatial lookup with the given positions and smoothing radius, and then you can query the grid cells around the current position instead of iterating over all the positions.

Here's an implementation of the CalculateDensity function using spatial grid optimization:

#include <cmath>

float SmoothingKernel(float smoothingRadius, float distance);

float CalculateDensity(const std::vector<std::pair<float, float>>& positions, float smoothingRadius, int index, const std::vector<Entry>& spatialLookup, const std::vector<int>& startIndices) {
	float density = 0;
	const float mass = 1;

	std::pair<float, float> position = positions[index];

	// Update the spatial lookup and start indices
	// UpdateSpatialLookup(positions, smoothingRadius, spatialLookup, startIndices);

	// Get the cell coordinates of the current position
	auto [cellX, cellY] = PositionToCellCoord(position, smoothingRadius);

	// Loop over neighboring cells
	for (int x = cellX - 1; x <= cellX + 1; ++x) {
		for (int y = cellY - 1; y <= cellY + 1; ++y) {

			// Get the start index for the current cell
			uint32_t cellKey = GetKeyFromHash(HashCell(x, y));
			int startIndex = startIndices[cellKey];

			// If the start index is not set, continue with the next cell
			if (startIndex == std::numeric_limits<int>::max()) {
				continue;
			}

			// Iterate over the points in the current cell
			for (int i = startIndex; i < spatialLookup.size() && spatialLookup[i].cellKey == cellKey; ++i) {
				std::pair<float, float> p = positions[spatialLookup[i].index];

				// Calculate the distance and influence
				float dx = p.first - position.first;
				float dy = p.second - position.second;
				float dst = std::sqrt(dx * dx + dy * dy);

				if (dst > smoothingRadius) continue;

				float influence = SmoothingKernel(smoothingRadius, dst);
				density += mass * influence;
			}
		}
	}

	return density;
}

Before calling the CalculateDensity function, make sure to call UpdateSpatialLookup with the positions and smoothing radius, and pass the resulting spatialLookup and startIndices vectors to the CalculateDensity function :

std::vector<Entry> spatialLookup;
std::vector<int> startIndices;
UpdateSpatialLookup(positions, smoothingRadius, spatialLookup, startIndices);

float density = CalculateDensity(positions, smoothingRadius, index, spatialLookup, startIndices);

This implementation should provide better performance by only iterating over the neighboring cells instead of all the positions.Make sure to implement the helper functions PositionToCellCoord, GetKeyFromHash, and HashCell if you haven't already.