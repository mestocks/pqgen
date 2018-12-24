#ifndef pq_htable_h__
#define pq_htable_h__

struct Hashll {
  char *key;
  void *value;
  struct Hashll *next;
};

struct HashTable {
  unsigned long size;
  struct Hashll **llptr;
};

extern unsigned long pq_FNV1a(char *key);
extern unsigned long pq_index_hash(char *key, unsigned long size);

extern void pq_create_hash(struct HashTable *hash, unsigned long size);
extern void pq_free_hash(struct HashTable *hash);
extern void pq_free_hash2(struct HashTable *hash);
extern void pq_insert_hash(struct HashTable *hash, char *key, void *value);
extern struct Hashll *pq_insert_hash2(struct HashTable *hash, char *key, void *value);
extern void *pq_lookup_hash(struct HashTable *hash, char *key);
extern void pq_print_hash(struct HashTable *hash);

extern void pq_multi_insert_hash(struct HashTable *hash, unsigned int argc, char **argv);

#endif
