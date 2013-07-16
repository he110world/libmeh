#ifndef COMPILER_H
#define COMPILER_H

// Regular expression
typedef struct NFAState_s{
	struct NFAState_s *next[2];
	int input[2];
	int i;
} NFAState_t;

typedef struct {
	NFAState_t *start, *final;
} NFA_t;


/* implement the Set as binary tree */
typedef struct BTree_s{
	struct BTree_s *less, *greater;
	int value[1];		/* variable sized */
} BTree_t;


typedef struct {
	BTree_t *btree;
	int extra;
} Set_t;


typedef struct listnode_s{
	struct listnode_s *prev, *next;
	int value;
} listnode_t;


typedef struct {
	listnode_t *head;
} List_t;


/* int key -> int value */
typedef struct Dict_s {
	BTree_t *btree;
} Dict_t;


typedef struct {
} Hash_t;

typedef struct {
	Dict_t **state;
	char *accept;
	int nstates;
} DFA_t;

typedef DFA_t RegEx_t;




#endif
