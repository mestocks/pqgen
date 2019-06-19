#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "pq_htable.h"

const uint64_t FNV_PRIME = 1099511628211;
const uint64_t FNV_OFFSET = 14695981039346656037;


unsigned long pq_FNV1a(char *key)
{  
  unsigned long i;
  unsigned long hash;
  hash = FNV_OFFSET;
  for (i = 0; i < strlen(key); i++) {
    hash ^= key[i];
    hash *= FNV_PRIME;
  }
  return hash;
}


unsigned long pq_index_hash(char *key, unsigned long size)
{
  unsigned long index;
  unsigned long hval;

  hval = pq_FNV1a(key);
  index = hval % size;
  return index;
}


void pq_create_hash(struct HashTable *hash, unsigned long hsize)
{
  hash->size = hsize;
  hash->llptr = calloc(hsize, sizeof *hash->llptr);
}

void pq_free_hash(struct HashTable *hash)
{
  int i;
  struct Hashll *curr;
  struct Hashll *prev;
  for (i = 0; i < hash->size; i++) {
    if (hash->llptr[i] != NULL) {
      curr = hash->llptr[i];
      while (curr->next != NULL) {
	prev = curr;
	curr = curr->next;
	free(prev->key);
	free(prev->value);
	free(prev);
      }
      free(curr->key);
      free(curr->value);
      free(curr);
    }
  }
  free(hash->llptr);
}

void pq_free_hash2(struct HashTable *hash)
{
  unsigned long i;
  struct Hashll *curr, *prev;
  
  for (i = 0; i < hash->size; i++) {
    if (hash->llptr[i] != NULL) {
      curr = hash->llptr[i];
      while (curr->next != NULL) {
	prev = curr;
	curr = curr->next;
	free(prev);
      }
      free(curr);
    }
  }
  free(hash->llptr);
}

/*
 * pq_insert_hash - insert key:value into hash table
 *
 * 1) Hash char for index
 * 2) if array[index] == NULL allocate linked list node
 * 3) else navigate linked list until node->next == NULL and allocate node 
 * 4) or node->key == key and assign value to node->value
 * 
 * Responsibility is on the caller to free memory (can be done through a call 
 * to pq_free_hash);
*/
void pq_insert_hash(struct HashTable *hash, char *key, void *value)
{
  unsigned long index;
  struct Hashll *curr;
  struct Hashll *newll;
  struct Hashll **pnext;
  
  index = pq_index_hash(key, hash->size);
  curr = hash->llptr[index];
  pnext = &hash->llptr[index];
  for (;;) {
    if (curr == NULL) {
	newll = malloc(sizeof (struct Hashll));
	newll->key = key;
	newll->value = value;
	newll->next = NULL;
	*pnext = newll;
	pnext = &newll->next;
	break;
    } else if (strcmp(curr->key, key) == 0) {
      free(curr->key);
      free(curr->value);
      curr->key = key;
      curr->value = value;
      break;
    } else {
      pnext = &curr->next;
      curr = curr->next;
    }
  }
}

struct Hashll *pq_insert_hash2(struct HashTable *hash, char *key, void *value)
{
  unsigned long index;
  struct Hashll *old, *new, *prev, *curr;

  new = calloc(1, sizeof *new);
  new->key = key;
  new->value = value;
  new->next = NULL;
  
  index = pq_index_hash(key, hash->size);
  curr = hash->llptr[index];

  if (curr == NULL) {
    old = curr;
    hash->llptr[index] = new;
  } else {

    while (1) {
      if (curr == NULL) {
	prev->next = new;
	old = curr;
	break;
      } else if (strcmp(curr->key, key) == 0) {
	new->next = curr->next;
	prev->next = new;
	old = curr;
	break;
      }
      prev = curr;
      curr = curr->next;
    }
  }
  
  return old;
}


/*
 * pq_lookup_hash - retrieve value from hash table with key
 *
 * Returns a NULL pointer if key is not found.
 */
void *pq_lookup_hash(struct HashTable *hash, char *key)
{
  unsigned int index;
  struct Hashll *curr;

  index = pq_index_hash(key, hash->size);
  curr = hash->llptr[index];
  while (curr != NULL && strcmp(curr->key, key) != 0) {
    curr = curr->next;
  }
  
  return curr->value;
}

void pq_print_hash(struct HashTable *hash)
{
  unsigned long i;
  struct Hashll *curr;
  for (i = 0; i < hash->size; i++) {
    curr = hash->llptr[i];
    printf("%lu", i);
    while (curr != NULL) {
      //printf(" %p", (void *)curr);
      printf(" %p", curr);
      printf("\ninner1 %p %d\n", curr, curr == NULL);
      curr = curr->next;
      printf("\ninner2 %p %d\n", curr, curr == NULL);
    }
    printf("\nouter");
    //    printf(" %p", curr);
    printf("\n");
  }
  
}


void pq_multi_insert_hash(struct HashTable *arghash, unsigned int argc, char **argv)
{
  char *key;
  char *value;
  unsigned int i;
  for (i = 0; i < argc; i+=2) {
    key = malloc(128 * sizeof(char));
    value = malloc(128 * sizeof(char));
    strcpy(key, argv[i]);
    strcpy(value, argv[i+1]);
    pq_insert_hash(arghash, key, value);
  }
}

