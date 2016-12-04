#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char const *argv[])
{
	vector<int> v1;
	vector<int> v2;
	v1.push_back(1);
	v1.push_back(2);
	v1.push_back(3);

	v2.push_back(4);
	v2.push_back(5);

	vector<int> v3;

	v3 = v1;
	v1 = v2;
	v2 = v3;

	for(int i = 0; i < v1.size(); i++){
		cout<<v1[i];
	}
	return 0;
}