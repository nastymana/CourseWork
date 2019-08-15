// pch.cpp: исходный файл, соответствующий предкомпилированному заголовку; нужен для компиляции

#include "pch.h"

std::ofstream logfile;
std::ifstream loadfile;






//////****************


int findElement(const vector<int> &a, int st, int end, int value)
{
	int newvalue = -1;
	// binary search
	for (int i = st; i < end; i++)
	{
		if (a[i] == value)
		{
			newvalue = i;
		}
	}
	return newvalue;
}



// В целом этот файл можно пропустить, но не удаляйте его, если вы используете предкомпилированные заголовки.
