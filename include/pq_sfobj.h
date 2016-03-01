#ifndef pq_sfobj_h__
#define pq_sfobj_h__


// Watterson's theta

struct pqThetaW {
  int nsam;
  long long int s;
  long long int nsites;
  void (*add)(struct pqThetaW *, int);
  double (*eval)(struct pqThetaW *);
  void (*reset)(struct pqThetaW *);
};

extern void pqThetaWAdd(struct pqThetaW *ptr, int s);
extern double pqThetaWEval(struct pqThetaW *ptr);
extern void pqThetaWInit(struct pqThetaW *ptr, int nsam);
extern void pqThetaWReset(struct pqThetaW *ptr);


// Tajima's theta

struct pqThetaPi {
  int nsam;
  long long int nsites;
  long long int pisum;
  void (*add)(struct pqThetaPi *, int);
  double (*eval)(struct pqThetaPi *);
  void (*reset)(struct pqThetaPi *);
};

extern void pqThetaPiAdd(struct pqThetaPi *ptr, int nminor);
extern double pqThetaPiEval(struct pqThetaPi *ptr);
extern void pqThetaPiInit(struct pqThetaPi *ptr, int nsam);
extern void pqThetaPiReset(struct pqThetaPi *ptr);

// Folded site frequency spectrum

/*

  pqNewFoldedSFS(int nsam)

 */

struct pqFoldedSFS {
  int nbins;
};

extern int pqSizeofFoldedSFS(int nsam);
//extern void rwkFoldedSFSAdd(struct rwkFoldedSFS *ptr, int nminor);

//extern rwkFoldedSFS(struct rwkFoldedSFS *ptr, int nsam);

#endif
