//
//  MyHash.h
//  
//
//  Created by Degui Zhi on 8/16/13.
//
//

#ifndef ____MyHash__
#define ____MyHash__

#include <iostream>
// #define TABLE_SIZE    100000000
using namespace std;
#include <string>




class HashEntryInt {
private:
    int key;
    int value;
public:
    HashEntryInt(int key, int value) {
        this->key = key;
        this->value = value;
    }
    
    int getKey() {
        return key;
    }
    
    int getValue() {
        return value;
    }
};

class HashMapInt {
private:
    int TABLE_SIZE;
    HashEntryInt **table;
public:
    int size;
    HashMapInt(int t) {
        TABLE_SIZE=t;
        table = new HashEntryInt*[t];
        size = 0;
        for (int i = 0; i < t; i++)
            table[i] = NULL;
    }
    
    int get(int key) {
        int hash = (key % TABLE_SIZE);
        while (table[hash] != NULL && table[hash]->getKey() != key)
            hash = (hash + 1) % TABLE_SIZE;
        if (table[hash] == NULL)
            return 0;
        else
            return table[hash]->getValue();
    }
    
    void put(int key, int value) {
        int hash = (key % TABLE_SIZE);
        while (table[hash] != NULL && table[hash]->getKey() != key)
            hash = (hash + 1) % TABLE_SIZE;
        if (table[hash] != NULL) {
            delete table[hash];
            size --;
        }
        table[hash] = new HashEntryInt(key, value);
        size ++;
    }
    
    ~HashMapInt() {
        for (int i = 0; i < TABLE_SIZE; i++)
            if (table[i] != NULL)
                delete table[i];
        delete[] table;
    }
};








class HashEntry
{
public:
    HashEntry(string Key1,  int Value1)
    {
        this->Key = Key1;
        this->Value = Value1;
    }
    string Key;
    int Value;
};

class HashMap {
private:
    HashEntry **table;
    int hashString(string Key);
    int table_size;
    int testmain();
public:
    int size;
     HashMap(int table_size);
    
    int get(string key);
    
    void put(string key, int value);
    
     ~HashMap();
};



#endif /* defined(____MyHash__) */
