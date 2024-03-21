#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

struct KeyValue {
    string key;
    string value;
};

void GetNameList(const char* fd, vector<KeyValue>& data) {
    ifstream file(fd);
    if (!file.is_open()) {
        cerr << "Error opening file" << endl;
        return;
    }

    string line;
    const string pattern = "initUcell";

    while (getline(file, line)) {
        KeyValue kv;

        if (line.substr(0, pattern.size()) == pattern) {
            line.erase(0, pattern.size());

            size_t pos = line.find_first_not_of(" \t");
            if (pos != string::npos) {
                line.erase(0, pos);
                pos = line.find_first_of(" \t");
                kv.key = line.substr(0, pos);
                line.erase(0, pos);
                kv.value = line;
                data.push_back(kv);
            }
        } else {
            size_t pos = line.find_first_of(" \t");
            if (pos != string::npos) {
                kv.key = line.substr(0, pos);
                line.erase(0, pos);
                kv.value = line;
                data.push_back(kv);
            }
        }
    }

    file.close();
}

int main() {
    vector<KeyValue> data;
    GetNameList("data.in", data);

    for (const auto& kv : data) {
        cout << kv.key << " " << kv.value << endl;
    }

    return 0;
}
