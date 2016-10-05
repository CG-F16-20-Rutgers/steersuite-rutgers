#include "obstacles/GJK_EPA.h"


SteerLib::GJK_EPA::GJK_EPA()
{
}

//--------------------------GJK---------------------------//

//GJK algorithm
bool SteerLib::GJK_EPA::GJK(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, std::vector<Util::Vector>& simplex)
{
	//initialize, pick random direction
	Util::Vector D(1, 0, 0);
	simplex.push_back(Support(_shapeA, _shapeB, D));
	Util::Vector reverse = -1 * D;
	//loop simplex
	for (;;)
	{
		simplex.push_back(Support(_shapeA, _shapeB, reverse));
		if (simplex.back() * reverse <= 0) {
			return false;
		}
		else {
			if (SteerLib::GJK_EPA::DoSimplex(reverse, simplex))
				return true;
		}
	}
}

//return index of furthest point of shape vertices along direction d (algorithm from class slides)
int SteerLib::GJK_EPA::furthestIndex(Util::Vector dir, const std::vector<Util::Vector>& _shapeA) {
	double highest = dir*_shapeA[0];
	int highestIndex = 0;

	for (unsigned int i = 1; i < _shapeA.size(); i++)
	{
		double dot = dir*_shapeA[i];
		if (dot>highest)
		{
			highest = dot;
			highestIndex = i;
		}
	}
	return highestIndex;
}

//support function to calc minkowski difference
Util::Vector SteerLib::GJK_EPA::Support(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, Util::Vector dir)
{
	Util::Vector aPoint = _shapeA[SteerLib::GJK_EPA::furthestIndex(dir, _shapeA)];
	Util::Vector bPoint = _shapeB[SteerLib::GJK_EPA::furthestIndex(-1 * dir, _shapeB)];
	Util::Vector MinkowskiDifference = aPoint - bPoint;

	return MinkowskiDifference;
}

//check if points contain origin
bool SteerLib::GJK_EPA::DoSimplex(Util::Vector& dir, std::vector<Util::Vector>& simplex)
{
	Util::Vector p1 = simplex.back();
	Util::Vector p1o = -1 * p1;
	Util::Vector p2; Util::Vector p3; Util::Vector p12; Util::Vector p13;
	//triangle
	if (simplex.size() == 3)
	{
		p2 = simplex[1];
		p3 = simplex[0];
		p12 = p2 - p1;
		p13 = p3 - p1;

		dir = Util::Vector(p12.z, p12.y, -1 * p12.x);
		if (dir * p3 > 0){
			dir = dir * -1;
		}

		if (dir * p1o > 0){
			simplex.erase(simplex.begin() + 0);
			return false;
		}

		dir = Util::Vector(p13.z, p13.y, -1 * p13.x);
		if (dir * p1o > 0){
			simplex.erase(simplex.begin() + 1);
			return false;
		}
		return true;
	}else{ //line segment
		p2 = simplex[0];
		p12 = p2 - p1;

		dir = Util::Vector(p12.z, p12.y, -1 * p12.x);
		if (dir * p1o < 0){
			dir = -1 * dir;
		}
	}
	return false;
}



//--------------------------EPA---------------------------//

bool SteerLib::GJK_EPA::EPA(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, std::vector<Util::Vector>& simplex, float& pdepth, Util::Vector& pvector)
{
	for(;;)
	{
		float dist;
		int index;
		Util::Vector norm;
		//get the edge
		getEdge(simplex, dist, norm, index);
		//then get the support in direction of the edge's normal
		Util::Vector support = SteerLib::GJK_EPA::Support(_shapeA, _shapeB, norm); 

		float d = support*norm;

		if (d - dist <= 0)
		{
			pvector = norm;
			pdepth = dist;
			return true;
		}else{
			simplex.insert(simplex.begin() + index, support);
		}
	}
}

//get edge for EPA
void SteerLib::GJK_EPA::getEdge(std::vector<Util::Vector>& simplex, float& d, Util::Vector& norm, int& index)
{
	d = FLT_MAX;

	for (int i = 0; i < simplex.size(); i++)
	{
		int j;
		if (i + 1 == simplex.size())
			j = 0;
		else
			j = i + 1;
		Util::Vector v1 = simplex[i];
		Util::Vector v2 = simplex[j];
		Util::Vector edge = v2 - v1;
		Util::Vector origin = v1;
		//triple product to get vector from edge towards origin
		Util::Vector n = origin*(edge*edge) - edge*(edge*origin); 
		//normalize
		n = n / sqrt(pow(n.x, 2) + pow(n.y, 2) + pow(n.z, 2)); //normalize
		//distance from origin to edge
		float dist = n*v1; 
		if (dist < d)
		{
			d = dist;
			index = j;
			norm = n;
		}
	}
}

//Look at the GJK_EPA.h header file for documentation and instructions
bool SteerLib::GJK_EPA::intersect(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB)
{
	std::vector<Util::Vector> simplex;

	if (GJK(_shapeA, _shapeB, simplex))
	{
		EPA(_shapeA, _shapeB, simplex, return_penetration_depth, return_penetration_vector);
		return true;
	}
	return false; // There is no collision
}