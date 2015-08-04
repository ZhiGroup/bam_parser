//
//  MyHash.cpp
//  
//
//  Created by Degui Zhi on 8/16/13.
//
//

#include "MyHash.h"




HashMap::HashMap(int new_table_size) {
    table_size = new_table_size;
    size = 0;
    table = new HashEntry*[table_size];
    for (int i = 0; i < table_size; i++)
        table[i] = NULL;
}

int HashMap::get(string key) {
    int hash = hashString(key);
    while (table[hash] != NULL && table[hash]->Key != key)
        hash = (hash + 1) % table_size;
    if (table[hash] == NULL)
        return 0;
    else
        return table[hash]->Value;
}

void HashMap::put(string key, int value) {
    int hash = hashString(key);
    while (table[hash] != NULL && table[hash]->Key != key)
        hash = (hash + 1) % table_size;
    if (table[hash] != NULL)
    {
        delete table[hash];
        size --;
    }
    table[hash] = new HashEntry(key, value);
    size++;
}

 HashMap::~HashMap() {
    for (int i = 0; i < table_size; i++)
        if (table[i] != NULL)
            delete table[i];
    delete[] table;
}


int HashMap::hashString(string Key)
{//start hashString
    long h = 0;
    for(int i=0; i<Key.length(); i++)
    {
        //To get almost fair distributions of nodes over the array
        h = (h << 3) ^ Key[i];
    }
    return abs(h % table_size );
}//end hashString


int HashMap::testmain()
{
    HashMap h(128);
    cout << "start" << endl;
    h.put("A",1);
    cout << "ok" << endl;
    h.put("C",2);
    h.put("G",3);
    h.put("T",4);

    for (int i=0;i<1000000;i++)
        h.put("A",h.get("A")+1);
    cout << "hash[" << "A" << "]=" << h.get("A") << endl;
    return 0;
}

