#include"TElement.h"

TDrift::TDrift(string str2) {
	sscanf(str2.c_str(), "%f%f%f", &l, &r, &ry);
}

TField_Map::TField_Map(string str2) {
	char temp[15];
	sscanf(str2.c_str(), "%d%f%f%f%f%f%f%f%s", &geom, &l, &pha, &r, &kb, &ke,
			&ki, &ka, temp);
	filename = temp;
}
