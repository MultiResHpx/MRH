#ifndef MEASUREMENT_H
#define MEASUREMENT_H

class Measurement
{
public:
	Measurement() {x=0.0; y=0.0;}
	Measurement(float a,float b) {x = a; y = b;}

	inline float getX();
	inline float getY();
	inline void  setX(float);
	inline void  setY(float);

private:
	float x;
	float y;
};

inline float Measurement::getX()
{
	return x;
}

inline float Measurement::getY()
{
	return y;
}

inline void Measurement::setX(float a)
{
	x = a;
}

inline void Measurement::setY(float b)
{
	y = b;
}

#endif