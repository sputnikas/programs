#include <iostream>
#include <cstdio>
#include <cstring>

class String {
public:
    String(const char *s = "");
    String(const String &s);
    ~String();
    //
    friend std::ostream &operator<<(std::ostream &f, const String &s);
    // 
    const char* c_str() const;
    String operator +(const String &s);
    String& operator =(const String &s);
    String operator +=(const String &s);
    bool operator ==(const String &s);
    bool operator !=(const String &s);
    String operator !();
    bool operator <(const String &s);
    bool operator >(const String &s);
    bool operator <=(const String &s);
    bool operator >=(const String &s);
    char &operator [](int i);
    char operator [](int i) const;
    String operator ()(int left, int slength = 0);
private:
    char *data;
    int length;
};

String::String(const char *s) : length( (s != 0) ? strlen(s) : 0 )
{
    data = new char[length + 1];
    if (s != 0) {
        strcpy(data, s);
    } else {
        data[0] = '\0';
    }
}

String::String(const String &s) : length(s.length)
{
    data = new char[length + 1];
    if (s.data != 0) {
        strcpy(data, s.data);
    } else {
        data[0] = '\0';
    }
}

String::~String()
{
    delete [] data;
}

const char* String::c_str() const
{
    return data;
}

String String::operator +(const String &s)
{
    int new_length = length + s.length;
    char *tmp = new char[new_length + 1];
    for (int i = 0; i<new_length; i++) {
        if (i < length) {
            tmp[i] = data[i];
        } else {
            tmp[i] = s.data[i - length];
        }
    }
    tmp[new_length] = '\0';
    return String(tmp);
}

String String::operator +=(const String &s)
{
    int new_length = length + s.length;
    char *tmp = new char[new_length + 1];
    for (int i = 0; i<new_length; i++) {
        if (i < length) {
            tmp[i] = data[i];
        } else {
            tmp[i] = s.data[i - length];
        }
    }
    tmp[new_length] = '\0';
    delete [] data;
    data = tmp;
    length = new_length;
    return *this;
}

String& String::operator =(const String &s)
{
    if (this == &s) {
        return *this;
    }
    if (length != s.length) {
        delete [] data;
        length = s.length;
        data = new char[length + 1];
    }
    for (int i = 0; i<length; i++) {
        data[i] = s.data[i];
    }
    data[length] = '\0';
    return *this;
}

bool String::operator ==(const String &s)
{
    if (length != s.length) {
        return false;
    }
    for (int i = 0; i<length; i++) {
        if (data[i] != s.data[i]) {
            return false;
        }
    }
    return true;
}

bool String::operator !=(const String &s)
{
    return !(*this == s);
}

String String::operator !()
{
    char *tmp = new char[length + 1];
    for (int i = 0; i<length; i++) {
        if (isupper(data[i])) {
            tmp[i] = tolower(data[i]);
        } else if (islower(data[i])) {
            tmp[i] = toupper(data[i]);
        } else {
            tmp[i] = data[i];
        }
    }
    tmp[length] = '\0';
    return String(tmp);
}

bool String::operator <(const String &s)
{
    int l = (length > s.length) ? s.length : length;
    for (int i = 0; i<l; i++) {
        if (data[i] < s.data[i]) {
            return true;
        } else if (data[i] > s.data[i]) {
            return false;
        }
    }
    return false;
}

bool String::operator >(const String &s)
{
    int l = (length > s.length) ? s.length : length;
    for (int i = 0; i<l; i++) {
        if (data[i] < s.data[i]) {
            return false;
        } else if (data[i] > s.data[i]) {
            return true;
        }
    }
    return false;
}

bool String::operator <=(const String &s)
{
    int l = (length > s.length) ? s.length : length;
    for (int i = 0; i<l; i++) {
        if (data[i] < s.data[i]) {
            return true;
        } else if (data[i] > s.data[i]) {
            return false;
        }
    }
    return (s.length >= length);
}

bool String::operator >=(const String &s)
{
    int l = (length > s.length) ? s.length : length;
    for (int i = 0; i<l; i++) {
        if (data[i] < s.data[i]) {
            return false;
        } else if (data[i] > s.data[i]) {
            return true;
        }
    }
    return (length >= s.length);
}

String String::operator ()(int left, int slength)
{
    char *tmp = new char[slength + 1];
    for (int i = 0; i<slength; i++) {
        if (left + i < 0) {
            tmp[i] = data[(left + i) % length + length];
        } else {
            tmp[i] = data[(left + i) % length];
        }
    }
    tmp[slength] = '\0';
    return String(tmp);
}

char String::operator [](int i) const
{
    if (i < 0) {
        return data[i % length + length];
    } else {
        return data[i % length];
    }
}

char& String::operator [](int i)
{
    if (i < 0) {
        return data[i % length + length];
    } else {
        return data[i % length];
    }
}

std::ostream &operator<<(std::ostream &f, const String &s) {
    f << "\'" << s.data << "\'";
    return f;
}

void testString() {
    String s1, s2, s3;
    s1 = "char string cast";
    s2 = String("main constructor");
    s3 = s2;
    std::cout << "s1 = " << s1 << std::endl;
    std::cout << "s2 = " << s2 << std::endl;
    std::cout << "s3 = " << s3 << std::endl;
    std::cout << "s3.c_str() = "<< s3.c_str() << std::endl;
    std::cout << "s1 + s2 = " << s1 + s2 << std::endl;
    String s4 = "1212";
    s4 += s1;
    std::cout << "1212 + s1 = "<< s4 << std::endl;
    std::cout << "s3 == s2 ? " << (s3 == s2) << std::endl;
    std::cout << "s3 != s2 ? " << (s3 != s2) << std::endl;
    std::cout << "!s1 =  " << !(s1) << std::endl;
    std::cout << "s4 > s2 ? " << (s4 > s2) << std::endl;
    std::cout << "s3 <= s2 ? " << (s3 <= s2) << std::endl;
    std::cout << "s3 >= s2 ? " << (s3 >= s2) << std::endl;
    std::cout << "s3 < s2 ? " << (s3 < s2) << std::endl;
    std::cout << "s4[2] = " << s4[2] << std::endl;
    s4[2] = 'r';
    std::cout << "s4 after s4[2] = r " << s4 << std::endl;
    std::cout << "s4(20, 50) = " << s4(20, 50) << std::endl;
    std::cout << "s4[-1] = " << s4[-1] << std::endl;
}

int main() {
    testString();
    return 0;
}