/* hahs-table.c 
 * 
 * Creates a dynamically allocated hash table that will be used by sim-outorder.
 * The sim will use the hashes to keep track of baddr and btarget values in the
 * FWD flow for use in the REV flow.
 *
 * Code originated from https://www.w3resource.com/c-programming-exercises/hash/c-hash-exercises-4.php
 * 
 */
 
// Including necessary header files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TABLE_SIZE 10
#define LOAD_FACTOR_THRESHOLD 0.7

// Forward declarations
struct HashTable* createHashTable(int size);
void resizeHashTable(struct HashTable* hashTable, int newSize);
void insert(struct HashTable* hashTable, const char* key, const char* value);
const char* retrieve(struct HashTable* hashTable, const char* key);
void displayHashTable(struct HashTable* hashTable);
void freeHashTable(struct HashTable* hashTable);

// Structure for a key-value pair
struct KeyValuePair {
    char key[50];
    char value[50];
    struct KeyValuePair* next;
};

// Structure for the hash table
struct HashTable {
    int size;
    int itemCount;
    struct KeyValuePair** table;
};
// Hash function for strings (djb2 algorithm)
unsigned long hashFunction(const char* str) {
    unsigned long hash = 5381;
    int c;

    while ((c = *str++) != '\0') {
        hash = ((hash << 5) + hash) + c; // hash * 33 + c
    }

    return hash;
}

// Function to create a new key-value pair
struct KeyValuePair* createKeyValuePair(const char* key, const char* value) {
    struct KeyValuePair* newPair = (struct KeyValuePair*)malloc(sizeof(struct KeyValuePair));
    if (newPair != NULL) {
        strcpy(newPair->key, key);
        strcpy(newPair->value, value);
        newPair->next = NULL;
    }
    return newPair;
}

// Function to create a new hash table
struct HashTable* createHashTable(int size) {
    struct HashTable* newTable = (struct HashTable*)malloc(sizeof(struct HashTable));
    if (newTable != NULL) {
        newTable->size = size;
        newTable->itemCount = 0;
        newTable->table = (struct KeyValuePair**)calloc(size, sizeof(struct KeyValuePair*));
    }
    return newTable;
}

// Function to insert a key-value pair into the hash table
void insert(struct HashTable* hashTable, const char* key, const char* value) {
    unsigned long index = hashFunction(key) % hashTable->size;

    // Create a new key-value pair
    struct KeyValuePair* newPair = createKeyValuePair(key, value);

    // Insert the new pair at the beginning of the linked list
    newPair->next = hashTable->table[index];
    hashTable->table[index] = newPair;

    // Update the item count
    hashTable->itemCount++;

    // Check if resizing is needed
    if ((double)hashTable->itemCount / hashTable->size > LOAD_FACTOR_THRESHOLD) {
        resizeHashTable(hashTable, hashTable->size * 2); // Resize to double the current size
    }
}

// Function to retrieve the value associated with a key
const char* retrieve(struct HashTable* hashTable, const char* key) {
    unsigned long index = hashFunction(key) % hashTable->size;
    struct KeyValuePair* current = hashTable->table[index];

    // Traverse the linked list at the index
    while (current != NULL) {
        if (strcmp(current->key, key) == 0) {
            return current->value; // Key found, return the value
        }
        current = current->next;
    }
    return "Key not found"; // Key not found
}

// Function to display the contents of the hash table
void displayHashTable(struct HashTable* hashTable) {
    for (int i = 0; i < hashTable->size; i++) {
        printf("[%d] -> ", i);

        struct KeyValuePair* current = hashTable->table[i];
        while (current != NULL) {
            printf("(%s, %s) -> ", current->key, current->value);
            current = current->next;
        }

        printf("NULL\n");
    }
}

// Function to free the memory allocated for the hash table
void freeHashTable(struct HashTable* hashTable) {
    for (int i = 0; i < hashTable->size; i++) {
        struct KeyValuePair* current = hashTable->table[i];
        while (current != NULL) {
            struct KeyValuePair* temp = current;
            current = current->next;
            free(temp);
        }
    }

    free(hashTable->table);
    free(hashTable);
}

// Function to resize the hash table
void resizeHashTable(struct HashTable* hashTable, int newSize) {
    // Create a new hash table with the specified size
    struct HashTable* newTable = createHashTable(newSize);

    // Rehash and insert all existing key-value pairs into the new table
    for (int i = 0; i < hashTable->size; i++) {
        struct KeyValuePair* current = hashTable->table[i];
        while (current != NULL) {
            insert(newTable, current->key, current->value);
            current = current->next;
        }
    }

    // Update the hash table with the new size and table
    hashTable->size = newSize;
    free(hashTable->table);
    hashTable->table = newTable->table;

    // Free the memory allocated for the new table structure
    free(newTable);
}
