#include "compiler.h"


// Regular expression
/* Regular Expression Parser */
/* copyrighted by Zxwyeah */
/*
 * feb,7,2007 - feb,15,2007: started learning regexp - implemented the parser
 *
 *
 * todo:
 * Error();
 * Preprocess();
 */

/*
  NFA:
  (start) --symbol--> (end)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "regexp.h"


static void RE2NFA();
static void NFA2DFA();


static NFA_t NFAStack[1024];
static int sp = 0;

static char vocabenc[128];	/* whether the input symbol is encountered. */
static char vocab[128];
static int numvocabs;

static void Push( char input ){
    NFAState_t *s0, *s1;
    NFA_t *n;

    printf( "push %c\n", input );

    s0 = malloc( sizeof( NFAState_t ) );
    s1 = calloc( 1, sizeof( NFAState_t ) );

    s0->next[0] = s1;
    s0->next[1] = NULL;
    s0->input[0] = input;
    s0->i = 0;

    n = NFAStack + sp++;
    n->start = s0;
    n->final = s1;
}


static void Concat( ){
    NFA_t *n0, *n1;

    printf( "concat\n" );

    n0 = NFAStack + sp - 2;
    n1 = NFAStack + sp - 1;

    n0->final->next[0] = n1->start;
    n0->final->input[0] = -1;
    n0->final = n1->final;

    sp--;
}


static void Union( ){
    NFA_t *n0, *n1;
    NFAState_t *s0, *s1;

    n0 = NFAStack + sp - 2;
    n1 = NFAStack + sp - 1;

    s0 = malloc( sizeof( NFAState_t ) );
    s1 = calloc( 1, sizeof( NFAState_t ) );

    s0->next[0] = n0->start;
    s0->next[1] = n1->start;
    s0->input[0] = s0->input[1] = -1;
    s0->i = 0;

    n0->final->next[0] = n1->final->next[0] = s1;
    n0->final->input[0] = n1->final->input[0] = -1;

    n0->start = s0;
    n0->final = s1;

    sp--;

}


/* * is the closure operator */
static void Star( ){
    NFAState_t *s0, *s1;
    NFA_t *n;

    printf( "star\n" );

    n = NFAStack + sp - 1;

    s0 = malloc( sizeof( NFAState_t ) );
    s1 = calloc( 1, sizeof( NFAState_t ) );

    s0->next[0] = n->start;
    s0->next[1] = s1;
    s0->input[0] = s0->input[1] = -1;
    s0->i = 0;

    n->final->next[0] = n->start;
    n->final->input[0] = -1;

    n->final->next[1] = s1;
    n->final->input[1] = -1;

    n->start = s0;
    n->final = s1;
}


static NFA_t Pop( ){
    sp--;

    if( sp < 0 )
	sp = 0;

    return NFAStack[sp==0 ? 0 : sp-1];
}


static NFA_t *Top( ){
    return NFAStack + sp - 1;
}

/*
  Grammar of Regular Expression

  precedence: closure > concatenation > union

  regExp ::= unionExp

  unionExp ::= concateExp {"|" concateExp}* 

  concateExp ::= closureExp |
  closureExp {closureExp}*

  closureExp ::= atom |
  atom "*"

  atom ::= char |
  "(" regExp ")"

*/

static char *src;
static char tok;
static char *cur;

static void GetToken( ){
    tok = *cur++;
    printf( "%c ", tok );
}

static void Error( ){
    printf( "parse error\n" );
    exit( 1 );
}

static void Expect( char c ){
    if( tok != c ){
	printf( "current char is %i. expecting %c\n", tok, c );
	Error( );
    }

    tok = *cur++;
}


static void UnionExp( );

static void Atom( ){
    if( !tok )
	return;

    if( tok == '(' ){
	GetToken( );
	UnionExp( );
	Expect( ')' );
    }
    else {
	if( tok == ')' )
	    return;

	Push( tok );
	if( !vocabenc[tok] ){
	    vocabenc[tok] = 1;
	    vocab[numvocabs++] = tok;
	}
	GetToken( );
    }
}


static void ClosureExp( ){
    Atom( );

    if( tok == '*' ){
	Star( );
	GetToken( );
    }
}


static void ConcateExp( ){
    ClosureExp( );

    while( tok && tok != '|' && tok != ')' ){
	ClosureExp( );
	Concat( );
    }
}


static void UnionExp( ){
    ConcateExp( );

    while( tok == '|' ){
	GetToken( );
	ConcateExp( );
	Union( );
    }
}


/* produces a NFA */
static void Parse( char *string ){
    char *c;

    src = strdup( string );
    cur = src;

    GetToken();
    UnionExp( );

    free( src );
}


/* n is the number of extra data in int */
Set_t *SET_New( ){
    Set_t *s;

    s = malloc( sizeof( Set_t ) );
    s->btree = NULL;

    return s;
}

#define newbtree( n ) malloc( &((( BTree_t * )0)->value[n]) )


static BTree_t *btree_dupr( BTree_t *t, int n ){
    BTree_t *d;

    if( !t )
	return NULL;

    d = newbtree( n );
    memcpy( d->value, t->value, n*sizeof( int ) );
    d->less = btree_dupr( t->less, n );
    d->greater = btree_dupr( t->greater, n );

    return d;
}


Set_t *SET_Dup( Set_t *s ){
    Set_t *ds;

    ds = SET_New( );
    ds->btree = btree_dupr( s->btree, 1 );
    /* ds->extra = s->extra */

    return ds;
}

static void btree_insertr( BTree_t *t, int value[], int n ){
    if( *value < *t->value ){
	if( !t->less ){
	    t->less = newbtree( n );
	    t->less->less = t->less->greater = NULL;
	    memcpy( t->less->value, value, n*sizeof( int ) );
	}
	else {
	    btree_insertr( t->less, value, n );
	}
    }
    else if( *value > *t->value ){
	if( !t->greater ){
	    t->greater = newbtree( n );
	    t->greater->less = t->greater->greater = NULL;
	    memcpy( t->greater->value, value, n*sizeof( int ) );
	}
	else {
	    btree_insertr( t->greater, value, n );
	}
    }
    else {			/* the value already exists */
	memcpy( t->value, value, n*sizeof( int ) );
    }

}



/*
 *
 *
 *
 *
 */
void SET_Insert( Set_t *s, int value ){
    if( !s->btree ){
	s->btree = newbtree( 1 );
	s->btree->less = s->btree->greater = NULL;
	*s->btree->value = value;
    }
    else {
	btree_insertr( s->btree, &value, 1 );
    }
}


/* Only compares value to btree->value[0] */
static BTree_t *btree_findr( BTree_t *t, int value ){
    if( !t )
	return NULL;
    if( value < *t->value )
	return btree_findr( t->less, value );
    else if( value > *t->value )
	return btree_findr( t->greater, value );
    else
	return t;

}


static BTree_t *btree_find( BTree_t *tree, int value, BTree_t **parent ){
    BTree_t *t, *oldt;

    t = tree;
    oldt = NULL;
    while( t ){
	if( value < *t->value ){
	    oldt = t;
	    t = t->less;
	}
	else if( value > *t->value ){
	    oldt = t;
	    t = t->greater;
	}
	else{
	    *parent = oldt;
	    return t;
	}
    }

    return NULL;
}


/*
  If a node has one child, then replace the node with its child
  If a node has two children, then replace the value of the node with the smallest value of the right tree.

*/
void SET_Delete( Set_t *s, int value ){
    BTree_t *p, *t, *subt, *mint, *oldmint;

    t = btree_find( s->btree, value, &p );

    if( t->less && t->greater ){
	mint = t->greater;
	oldmint = t;
	while( mint->less ){
	    oldmint = mint;
	    mint = mint->less;
	}
	*t->value = *mint->value;

	/* delete the smallest node in the right subtree */
	oldmint->less = mint->greater;
	free( mint );

    }
    else {
	subt = t->less ? t->less : t->greater;
	if( p ){
	    if( value < p->value )
		p->less = subt;
	    else
		p->greater = subt;
	}
	else {
	    s->btree = subt;
	}

	free( t );
			 
    }
}


/* compare whether two sets are equal */
static int btree_equalr( BTree_t *t0, BTree_t *t1 ){
    if( t0 && t1 ) {
	if( *t0->value == *t1->value ){
	    return btree_equalr( t0->less, t1->less ) && 
		btree_equalr( t0->greater, t1->greater );
	}
	else 
	    return 0;
    }
    else if( !t0 && !t1 ){
	return 1;
    }
    else
	return 0;
}


/*
 *
 *
 *
 *
 *
 */
int SET_Equal( Set_t *s0, Set_t *s1 ){
    int i;	
    i = btree_equalr( s0->btree, s1->btree );

    return i;
}


static void unionr( Set_t *s, BTree_t *t ){
    if( !t )
	return;

    SET_Insert( s, *t->value );
    unionr( s, t->less );
    unionr( s, t->greater );
	
}


/*
 *
 *
 *
 *
 */
void SET_Union( Set_t *s0, Set_t *s1 ){
    unionr( s0, s1->btree );
}


static void btree_destroyr( BTree_t *t ){
    if( !t )
	return;

    btree_destroyr( t->less );
    btree_destroyr( t->greater );
    free( t );
}


void SET_Destroy( Set_t *s ){
    btree_destroyr( s->btree );
}



static int nvisits;

static void eclosurer( NFAState_t *s, Set_t *states ){
    if( !s )
	return;

    if( s->i >= nvisits )
	return;

    SET_Insert( states, s );

    s->i++;

    printf( "s->input[0] = %i, s->input[1] = %i\n", s->input[0], s->input[1] );

    if( s->input[0] == -1 ){
	printf( "haha\n" );
	eclosurer( s->next[0], states );
    }
    if( s->input[1] == -1 ){
	printf( "hahaha\n" );
	eclosurer( s->next[1], states );
    }

}


static void btree_eclosurer( BTree_t *t, Set_t *states ){
    Set_t *s;

    if( !t )
	return;

    s = SET_New( );
    eclosurer( *t->value, s );

    SET_Union( states, s );
    SET_Destroy( s );

    btree_eclosurer( t->less, states );
    btree_eclosurer( t->greater, states );

}


static Set_t *EpsilonClosure( Set_t *s ){
    Set_t *states;
    BTree_t *t;
    int i;

    nvisits++;

    states = SET_New( );
    btree_eclosurer( s->btree, states );

    return states;
}


static void mover( BTree_t *t, char input, Set_t *reach ){
    NFAState_t *s;

    if( !t )
	return;

    s = *t->value;
    if( input == s->input[0] )
	SET_Insert( reach, s->next[0] );
    if( input == s->input[1] )
	SET_Insert( reach, s->next[1] );

    mover( t->less, input, reach );
    mover( t->greater, input, reach );
}


static Set_t *Move( Set_t *s, char input ){
    Set_t *reach;

    reach = SET_New( );

    mover( s->btree, input, reach );

    return reach;
}



List_t *LIST_New( ){
    return calloc( 1, sizeof( List_t ) );
}


void LIST_Insert( List_t *list, int value ){
    listnode_t *head;

    head = malloc( sizeof( listnode_t ) );
    head->value = value;
    head->prev = NULL;
    head->next = list->head;
    if( list->head )
	list->head->prev = head;
    list->head = head;
}


/* Move node n from list s(ource) to list d(estiny) */
void LIST_Move( listnode_t *n, List_t *s, List_t *d ){
    if( !n->prev )
	s->head = n->next;
    else
	n->prev->next = n->next;

    if( n->next )
	n->next->prev = n->prev;

    n->next = d->head;
    n->prev = NULL;

    if( n->next )
	n->next->prev = n;

    d->head = n;

}


void LIST_Destroy( List_t *list ){
    listnode_t *l, *next;

    l = list->head;
    while( l ){
	next = l->next;
	free( l );
	l = next;
    }

    free( list );
}


Dict_t *DICT_New( ){
    Dict_t *d;

    d = malloc( sizeof( Dict_t ) );
    d->btree = NULL;

    return d;
}



void DICT_Set( Dict_t *d, int key, int value ){
    int keyvalue[2];

    if( !d->btree ){
	d->btree = newbtree( 2 );
	d->btree->less = d->btree->greater = NULL;
	d->btree->value[0] = key;
	d->btree->value[1] = value;

	return;
    }

    keyvalue[0] = key;
    keyvalue[1] = value;

    btree_insertr( d->btree, keyvalue, 2 );
}

int DICT_Get( Dict_t *d, int key, int *value ){
    BTree_t *t;

    t =  btree_findr( d->btree, key );

    if( t ) {
	*value = t->value[1];
	return 1;
    }
    else {
	return 0;
    }

}

void DICT_Remove( Dict_t *d, int key ){
    SET_Delete( d, key );
}


void DICT_Destroy( Dict_t *d ){
    SET_Destroy( d );
}


typedef struct {
    Set_t *start;
    char input;
    Set_t *final;
} transition_t;


static int isacceptingr( BTree_t *t ){
    if( !t )
	return 0;

    if( !(( NFAState_t *)*t->value)->next[0] )
	return 1;

    return isacceptingr( t->less ) || isacceptingr( t->greater );
}

static int IsAccepting( Set_t *s ){
    return isacceptingr( s->btree );
}


static void deletenfar( NFAState_t *s, int n ){
    static int i = 0;

    if( !s )
	return;

    if( s->i == n )
	return;

    s->i = n;

    deletenfar( s->next[0], n );
    deletenfar( s->next[1], n );

    free( s );
}


static void DeleteNFA( NFAState_t * s ){
    nvisits++;

    deletenfar( s, nvisits );
}

/* convert NFA to DFA */
static DFA_t *SubsetConstruction( NFAState_t *s ){
    List_t *moved, *unmoved;
    Set_t *start;
    Set_t *mset, *umset;
    int i;
    listnode_t *um, *next, *m, *tr;
    int found;
    List_t *transition;
    transition_t *trans;
    int numstates;
    DFA_t *d;

    unmoved = LIST_New( );
    moved = LIST_New( );
    transition = LIST_New( );

    start = SET_New( );

    SET_Insert( start, s );
    printf( "hello = %i\n", (( NFAState_t * )start->btree->value[0])->input[0] );
    printf( "hello2 = %i\n", s->input[0] );

    LIST_Insert( unmoved, EpsilonClosure( start ) );

    numstates = 0;
    while( unmoved->head ){
	umset = unmoved->head->value;
	umset->extra = numstates; /* id of the DFA state */

	LIST_Move( unmoved->head, unmoved, moved );
	numstates++;

	for( i=0; i<numvocabs; i++ ){
	    mset = EpsilonClosure( Move( umset, vocab[i] ) );

	    printf( "mset = %i\n", mset );

	    if( !mset->btree ){
		free( mset );
		continue;
	    }

	    /* Do we already have this DFA state? */
	    found = 0;
	    for( m=moved->head; m; m=m->next ){
		if( SET_Equal( mset, m->value ) ){
		    found = 1;
		    SET_Destroy( mset );
		    mset = m->value;
		    break;
		}
	    }

	    /* If not, let's have it. */
	    if( !found ){
		LIST_Insert( unmoved, mset );
	    }

	    /* New state transition */
	    trans = malloc( sizeof( transition_t ) );
	    trans->start = umset;
	    trans->input = vocab[i];
	    trans->final = mset;
	    LIST_Insert( transition, trans );

	}
    }

    /* Now we've got a list of DFA states (sets of NFA states); A list of DFA transitions */
    /* Convert the DFAs to the final DFA transition table */
    d = malloc( sizeof( DFA_t ) );
    d->nstates = numstates;
    d->state = calloc( numstates, sizeof( Dict_t ) );
    d->accept = malloc( numstates );

    printf( "numstates = %i\n", numstates );
    for( i=0; i<numstates; i++ ){
	d->state[i] = DICT_New( );
    }

    m = moved->head;
    i = numstates - 1;
    while( m ){
	d->accept[i] = IsAccepting( m->value );
	printf( "accept = %i\n", d->accept[i] );
	printf( "extra = %i\n", ((Set_t *)m->value)->extra );
	m = m->next;
	i--;
    }

    tr = transition->head;
    printf( "transitions:\n" );
    while( tr ){
	trans = tr->value;

	printf( "%i, %c, %i\n", trans->start->extra, trans->input, trans->final->extra );
	DICT_Set( d->state[ trans->start->extra ], 
		  trans->input,
		  trans->final->extra );

	tr = tr->next;
    }

    /* delete all used stuff */

    LIST_Destroy( unmoved );

    /* DFA states */
    for( m=moved->head; m; m=next ){
	next = m->next;
	SET_Destroy( m->value );
    }

    LIST_Destroy( moved );

    /* NFA states */
    DeleteNFA( s );

    /* transition list */
    tr = transition->head;
    while( tr ){
	next = tr->next;
	free( tr->value );
	free( tr );
	tr = next;
    }

    return d;
}


/*
  Extended regular expression grammar
*/
void Preprocess( char *string ){

}


RegEx_t *RE_Compile( char *string ){
    memset( vocabenc, 0, 128 );
    nvisits = 0;
    Parse( string );

    printf( "sp = %i\n", Top()->start->input[0] );
    return SubsetConstruction( Top()->start );
}

int RE_Match( char *string, RegEx_t *re ){
    char *c;
    int s;

    s = 0;
    c = string;
    while( *c ){
	printf( "%i ", s );
	if( !DICT_Get( re->state[s], *c, &s ) )
	    return 0;

	c++;
    }

    return re->accept[s];
}


void RE_Delete( RegEx_t *re ){
    int i;

    for( i=0; i<re->nstates; i++ )
	DICT_Destroy( re->state[i] );

    free( re->accept );

    free( re );
}


#if 0
int main( ){
    RegEx_t *re;

    re = RE_Compile( "a|s(d)*a|f" );

    printf( "match = %i\n", RE_Match( "sdddd", re ) );
	
    RE_Delete( re );

    return 0;
}
#endif



// Compiler (Minic programming language)
/*
  minic:
  =============================
  LEXON:

  if
  else
  while
  for
  (
  )
  {
  }
  +
  -
  *
  /
  %
  =
  >
  <
  >=
  <=
  ==
  !=
  int
  float
  char * ==> string
  ,
  ;
  break
  continue
  return
  "

  =============================
  GRAMMAR:

  <translate-unit> = <extern-decl> { <extern-decl> }

  <extern-decl> = <type-spec> <identifier> <decl-tail>

  <type-spec> = "int" | "float" | "char" "*"

  <decl-tail> = 	"(" {<var-list>} ")" <func-decl-tail> |
  <var-decl-tail> ";"  |
  "=" <const> <var-decl-tail> ";"

  <var-decl-tail> = {"," <identifier> {"=" <const>} }

  <var-list> = <type-spec> <identifier> {"," <var-list> }

  <func-decl-tail> = "{" {<local-var-decl>} { <statement> } "}"

  <local-var-decl> = <type-spec> <identifier> ["=" <const>] <local-var-decl-tail> | $

  <local-var-decl-tail> = {"," <identifier> ["=" <const>] } ";"

  <statement> = 	<expr-statement>  |
  <compound-statement>  |
  <selection-statement>  |
  <iteration-statement>  |
  <jump-statement>

  <expr-statement> = <expr> ";"

  <compound-statement> = "{" {<statement>} "}"

  <selection-statement> = "if" "(" <expr> ")" <statement> <if-tail>

  <if-tail> = "else" <statement>  |  $

  <iteration-statement> = 	"while" "(" <expr> ")" <statement>

  <jump-statement> =		"continue"  ";"  |
  "break" ";"  |
  "return" [<expr] ";"

  <expr> ::= <assign> (0)
  <assign> ::= <l-or> { "=" <l-or> } (1)
  <l-or> ::= <l-and> { "||" <l-and> } (2)
  <l-and> ::= <equal> { "&&" <equal> } (3)
  <equal> ::= <rel> { ["=="|"!="] <rel> } (4)
  <rel> ::= <add> { ["<"|">"|"<="|">="] <add> } (5)
  <add> ::= <mul> { ["+"|"-"] <mul> } (6)
  <mul> ::= <unary> { ["*"|"/"|"*"] <unary> } (7)
  <unary> ::= {"+"|"-"|"!"} <primary> (8)
  <primary> ::= <id> { "(" <expr> ")" } (9)
  | <const>
  | "(" <expr> ")"

  | "[" <expr> "]"
  | <id> "." <id> { "(" <expr> ")" }



  <expr-list> = 	<expr> { "," <expr> }

  <constant> =	"int-const"  |
  "float-const" |
  "string-const"

  ====================================
  OPERATOR PRECEDENCE
  ()
  + - !
  * / %
  + -
  < <= > >=
  == !=
  &&
  ||
  =

  a = b || c && d == e < f + g * +h 

  a = b || c || d <=> (b || c) || d
  a = b || c && d <=> b || (c && d)
  a = b && c || d && e
  a + b
  a = 33
  a = b

  a + 33
  a + b
  33 + a
  33 + 33

  <assign> ::= <id> {"=" <assign> | <smth-else> <logic-or-expr>}
  | <logic-or-expr>	/* begin with "+" | "-" | "!" | <const> | "("

  <logic-or> ::= <logic-and> {"||" <logic-or> }



  ====================================
  SEMANTICS:

  variable name -> address
  statement -> operations
  function -> operations


  ====================================
  VIRTUAL MACHINE:

  program counter [pc], operand stack top [optop], frame pointer [fp], and variables pointer [vars]
  registers: pc, sp, fp

  _______________________________ STACK FRAME _______________________________________________________
  |                                                                                                   |
  \/__________________________________________________________________________________________________\/
  | ARGS...| CONTEXT (PC, SP, FP) | LOCAL VARS | OPERAND STACK (args of callee is on the stack, too)  |
  |________|______________________|____________|______________________________________________________|
  /\                                                                  /\  
  |                                                                   |
  fp                                                                  sp
  opcodes:

  loadpc, storepc,
  loadsp, storesp,
  loadfp, storefp,
  loadvar, storevar,

  call, ret,
  iloads, istores,
  iload, istore,
  imul, idiv, iadd, isub, imod,
  ig, il, ige, ile,
  ie, ine,  

  floads, fstores,
  fload, fstore,
  fmul, fdiv, fadd, fsub, fmod,
  fg, fl, fge, fle,
  fe, fne,

  jz, jnz,
  jmp,

  call,	// call native function


*/

#include <stdio.h>
#include <setjmp.h>
#include <stdlib.h>
#include <string.h>

#include "../include/compiler.h"

static FILE *f;

enum {TOK_IF = 256, TOK_ELSE, TOK_WHILE, TOK_BREAK, TOK_CONTINUE, TOK_RETURN, TOK_INT,
      TOK_FLOAT, TOK_FOR, TOK_VOID, TOK_GE, TOK_LE, TOK_EQ,
      TOK_NEQ, TOK_INT_CONST, TOK_FLOAT_CONST, TOK_STRING_CONST, TOK_END, TOK_ID,
      TOK_L_OR, TOK_L_AND};

enum {OP_LOADC, OP_LOAD, OP_STORE, OP_LOADG, OP_STOREG, OP_IADD, OP_ISUB,
      OP_IMUL, OP_IDIV, OP_MOD, OP_IGT, OP_ILT, OP_IGE, OP_ILE, OP_IEQ,
      OP_INE, OP_I2F, OP_F2I, OP_FADD, OP_FSUB, OP_FMUL, OP_FDIV, OP_FGT,
      OP_FLT, OP_FGE, OP_FLE, OP_FEQ, OP_FNE, OP_JZ, OP_JNZ, OP_CALL, 
      OP_RET, OP_RETN, OP_SWAP, OP_DUP, OP_COMPRESS, OP_POP, OP_L_OR, 
      OP_L_AND, OP_NEG, OP_END, OP_NOT, OP_JMP, OP_NOP, OP_CALLN, OP_ENLARGE};

static char *kwd[] = { "if", (char *)TOK_IF, 
		       "else", (char *)TOK_ELSE,
		       "while", (char *)TOK_WHILE,
		       "break", (char *)TOK_BREAK,
		       "continue", (char *)TOK_CONTINUE,
		       "return", (char *)TOK_RETURN,
		       "int", (char *)TOK_INT,
		       "float", (char *)TOK_FLOAT,
		       "for", (char *)TOK_FOR,
		       "void", (char *)TOK_VOID };

static int numkwd = sizeof (kwd) / 2 / sizeof (char *);
union {
    int i;
    float f;
    char *str;			/* the content of sem.str and sem.id will be overwritten after next gettok ()*/
    char *id;
} sem;

static char buf[1000];
static char *pbuf;
static int c;
static int tok;

/*
 */
static void getch () {
    c = fgetc (f);
    printf ("%c", c);
}
  
/* 
 */
jmp_buf jb;
int god;
static void prerror (char *err) {
    printf (" <- error: %s\n", err);
    longjmp (jb, 1);
}
  
/* 
 */
/************************************/
/* static int isnum (char ch) {     */
/*   return ch >= '0' && ch <= '9'; */
/* }                                */
/************************************/

#define isnum(ch) (ch >= '0' && ch <= '9')
/* 
 */
int gettok () {
    int numdots;
    int i;
    char c0;
  
    while (c == ' ' || c == '\n' || c == '\t' || c == '\r' || c == '/') {
	while (c == ' ' || c == '\n' || c == '\r' || c == '\t') {
	    getch ();
	}    
	while (c == '/') {
	    getch ();
	    if (c == '*') {
		getch ();
		while (1) {
		    if (c == EOF) {
			prerror ("unterminated comment");
		    }        
		    if (c == '*') {
			getch ();
			if (c == '/') {
			    getch ();
			    break;
			}
		    }
		    else {
			getch ();
		    }
		}
	    }
	    else {
		tok = '/';
		return tok;
	    }
	}    
    }

    if (c == EOF) {
	tok = TOK_END;
    }
    else {
	pbuf = buf;
	numdots = 0;
	if (c == '.' || isnum(c) ) { /* numbers */
	    do {
		*pbuf++ = c;
		if (c== '.') {
		    numdots++;
		}
		getch ();
	    } while (c == '.' || isnum(c) );
	    *pbuf = 0;
	    if (numdots == 0) {
		tok = TOK_INT_CONST;
		sscanf (buf, "%d", &sem.i);
		god++;
	    }
	    else if (numdots == 1) {
		tok = TOK_FLOAT_CONST;
		sscanf (buf, "%f", &sem.f);
	    }
	    else {
		prerror ("number has more than one dots");
	    }
	}
	else if (c == '_' || isalpha (c)) {	/* identifier or keyword */
	    do {
		*pbuf++ = c;
		getch ();
	    } while (isalnum (c) || c == '_');
	    *pbuf = 0;
	    for (i=0; i<2*numkwd; i+=2) {
		if (!strcmp (buf, kwd[i])) {
		    tok = (int) kwd[i+1];
		    break;
		}
	    }
	    if (i == 2*numkwd) {
		tok = TOK_ID;
		sem.id = buf;
	    }
	}
	else if (c == '\"') {
	    do {
		getch ();
		if (c == '\\') {
		    getch ();
		    switch (c) {
		    case 't':
			*pbuf++ = '\t';
			break;
		    case 'n':
			*pbuf++ = '\n';
			break;
		    case '\\':
			*pbuf++ = '\\';
			break;
		    case '\"':
			*pbuf++ = '\"';
			break;
		    default:
			prerror ("invalid string");
		    }
		    getch ();
		    continue;
		}
		*pbuf++ = c;
	    } while (c != '\"');
	    getch ();
	    *pbuf = 0;
	    tok = TOK_STRING_CONST;
	    sem.str = buf;
	}
	else {
	    switch (c) {
	    case '(':
	    case ')':
	    case '{':
	    case '}':
	    case '+':
	    case '-':
	    case '*':
	    case '%':
	    case ';':
	    case ',':
		tok = c;
		getch ();
		break;
	    case '>':
	    case '<':
	    case '=':
	    case '!':
		c0 = c;
		getch ();
		if (c == '=') {
		    switch (c0) {
		    case '>':
			tok = TOK_GE;
			break;
		    case '<':
			tok = TOK_LE;
			break;
		    case '=':
			tok = TOK_EQ;
			break;
		    case '!':
			tok = TOK_NEQ;
			break;
		    }
		    getch ();
		}
		else {
		    tok = c0;
		}
		break;
	    default:
		prerror ("unknown token");
		break;
	    }
      
	}
    }  

    return tok;
}

/* 
 */
char *gsym[10000];
int numgsyms;

char *lsym[100];
int numlsyms;

int gvar[10000];
int global = gvar;

int mem[100000];

int bytecode[100000];
int numops;

int pc, sp, fp;

int label[10000];
int numlabels;

char strmem[10000];
int numchars;

typedef struct {
    int islvalue;
    int gvar_addr;
    int lvar_addr;
    int isgvar;
    int type;
    int isfunc;
    int isdecl;
    int isbuildin;
    int numargs;
    int numlvars;
    int argtype[12];
} context_t;

context_t gvarcont[10000];
context_t lvarcont[100];


/* 
 */
static int genlab () {
    return numlabels++;
}

static void patchlab (int lab, int addr) {
    label[lab] = addr;
}

static void backpatchlab (int lab, int addr) {
    bytecode[label[lab]] = addr;
}

/* 
 */
static int addgsym (char *sym) {
    gsym[numgsyms++] = strdup (sym);
  
    return numgsyms-1;
}

static int addlsym (char *sym) {
    lsym[numlsyms++] = strdup (sym);
  
    return numlsyms-1;
}

/* 
 */
static int findgsym (char *sym) {
    int i;

    for (i=0; i<numgsyms; i++) {
	if (!strcmp (sym, gsym[i])) {
	    return i;
	}
    }
    if (i == numgsyms) {
	return -1;
    }
}

static int findlsym (char *sym) {
    int i;

    for (i=0; i<numlsyms; i++) {
	if (!strcmp (sym, lsym[i])) {
	    return i;
	}
    }
    return -1;
}

/*
  DELETE ME! (build in function test) 

  sptr -- stack pointer
*/
int print_api (int sptr) {
    printf ((char *)&(strmem[mem[sptr]]));

    return 0;
}


int printint_api(int sptr)
{
    printf("%i\n",mem[sptr]);
}


int printfloat_api(int sptr)
{
    printf("%f\n",*((float *)&(mem[sptr])));
}


/* 
 */
typedef void (* buildinfpvoid_t) (int);
typedef int (* buildinfp_t) (int);


typedef struct {
    char *name;
    buildinfp_t fp;
} buildin_t;

buildin_t buildin[1000] = {{"print", print_api}, {"printint", printint_api}, {"printfloat", printfloat_api}};
int numbuildins = 3;

static int findbuildin (char *func) {
    int i;

    for (i=0; i<numbuildins; i++) {
	if (!strcmp (buildin[i].name, func)) {
	    return i;
	}
    }

    return -1;
}

/* 
 */
static void resetlsym () {
    numlsyms = 0;
}


/* 
 */
static void emit (int b) {
    bytecode[numops++] = b;
}

/* 
 */
static void expect (int t) {
    if (tok != t) {
	prerror ("expect");
    }
    gettok ();
}

/* 
 */
static int addstr (char *str) {
    int old;

    old = numchars;
    numchars += strlen (str) + 1;
    strcpy (strmem+old, str);

    return old;
}


/* 
 */
typedef struct {
    int funcaddr;
    int opaddr;
} funcpatch_t;

funcpatch_t fpatch[10000];
int numfpatchs;
static void funcpatch (int funcaddr, int opaddr) {
    fpatch[numfpatchs].funcaddr = funcaddr;
    fpatch[numfpatchs++].opaddr = opaddr;
}

static void funcbackpatch () {
    int i;
    int fa, oa;

    for (i=0; i<numfpatchs; i++) {
	fa = fpatch[i].funcaddr;
	oa = fpatch[i].opaddr;
	if (gvarcont[fa].isdecl) {
	    prerror ("undefined reference to function");
	}
	bytecode [oa] = gvar [fa];
    }
}


/* 
   declearation
*/
static void block (int isfuncdecl, int funcaddr, int breakable);

char buf2[1024];
enum {TYPE_INT, TYPE_FLOAT, TYPE_VOID};

static void decl (int isglobal) {
    int type;
    int g, l;
    int straddr;
    int bldin;
    int i;
    int needdef;

    while (tok == TOK_INT || tok == TOK_FLOAT || tok == TOK_VOID) {
	switch (tok) {
	case TOK_INT:
	    type = TYPE_INT;
	    break;
	case TOK_FLOAT:
	    type = TYPE_FLOAT;
	    break;
	case TOK_VOID:
	    type = TYPE_VOID;
	    break;
	}
	gettok ();
	while (tok != ';') {
	    needdef = 0;
	    if (tok == ',') {
		gettok ();
	    }
	    if (tok != TOK_ID) {
		prerror ("invalid decl");
	    }
	    if (isglobal) {
		if ((g = findgsym (sem.id)) != -1) {
		    if (gvarcont[g].isdecl) {
			needdef = 1;
		    }
		    else {
			prerror ("global var redefine");
		    }
		}        
		else {
		    bldin = findbuildin (sem.id);
		    g = addgsym (sem.id);
		    gvarcont[g].type = type;
		}
	    }
	    else {
		if (findlsym (sem.id) != -1) {
		    prerror ("local var redefine");
		}
		l = addlsym (sem.id);
		lvarcont[l].type = type;
	    }
	    gettok ();		/* id */
	    if (tok == '(') {
		if (!isglobal) {
		    prerror ("local function decl is not support");
		}
		resetlsym ();
		gettok ();		/* ( */
		gvarcont[g].isfunc = 1;
		gvar[g] = numops;
		i = 0;
		while (tok != ')') {
		    if (tok == ',') {
			gettok ();
		    }
		    if ((tok != TOK_INT) && (tok != TOK_FLOAT)) {
			prerror ("error argument type");
		    }
		    switch (tok) {
		    case TOK_INT:
			gvarcont[g].argtype[i] = TYPE_INT;
			break;
		    case TOK_FLOAT:
			gvarcont[g].argtype[i] = TYPE_FLOAT;
			break;
		    }
		    gettok ();		/* type */
		    if (tok != TOK_ID) {
			prerror ("argument has no name");
		    }
		    if (findlsym (sem.id) != -1) {
			prerror ("argument redefine");
		    }
		    addlsym (sem.id);
		    i++;
		    gettok ();		/* id */
		}
		gettok ();		/* ) */
		if (tok == '{') {	/* function defination */
		    if (bldin != -1) {
			prerror ("conflict with build-in function");
		    }	    
		    gvarcont[g].isdecl = 0;
		    gvarcont[g].numargs = i;
		    gvarcont[g].isbuildin = 0;
		    /* Be careful! Addressing of local var is different from global var */
		    for (i=0; i<gvarcont[g].numargs; i++) {
			lvarcont[i].lvar_addr = -(1 + gvarcont[g].numargs - i); /* addressing part 1: argument addr */
			lvarcont[i].type=gvarcont[g].argtype[i];
		    }        
		    block (1, g, 0);
		    break;
		}
		else {			/* function declaration */
		    if (needdef) {
			prerror ("function protype redefine");
		    }
		    if (bldin != -1) {
			gvarcont[g].isdecl = 0;
			gvarcont[g].isbuildin = 1;
			gvar[g] = buildin[g].fp;
			gvarcont[g].numargs = i;
		    }
		    else {
			gvarcont[g].isdecl = 1;
			gvarcont[g].isbuildin = 0;
			gvarcont[g].numargs = i;
		    }
		}
	    }
	    else if (tok == '=') {
		if (type == TYPE_VOID) {
		    prerror ("storage size isn't known");
		}
		gettok ();		/* = */
		if (type == TYPE_INT) {
		    if (tok == TOK_INT_CONST) {
			if (isglobal) {
			    gvar[g] = sem.i;
			}
			else {
			    emit (OP_LOADC);
			    emit (sem.i);
			    emit (OP_STORE);
			    emit (l+1);
			}
		    }
		    else if (tok == TOK_FLOAT_CONST) {
			if (isglobal) {
			    gvar[g] = sem.f;
			}
			else {
			    emit (OP_LOADC);
			    emit (sem.i);
			    emit (OP_F2I);
			    emit (OP_STORE);
			    emit (l+1);
			}
		    }
		    else if (tok == TOK_STRING_CONST) {
			straddr = addstr (sem.str);
			if (isglobal) {
			    gvar[g] = straddr;
			}
			else {
			    emit (OP_LOADC);
			    emit (straddr);
			    emit (OP_STORE);
			    emit (l+1);
			}
		    }
		}
		else if (type == TYPE_FLOAT) {
		    if (tok == TOK_FLOAT_CONST) {
			if (isglobal) {
			    gvar[g] = sem.i;
			}
			else {
			    emit (OP_LOADC);
			    emit (sem.i);
			    emit (OP_STORE);
			    emit (l+1);
			}
		    }
		    else if (tok == TOK_INT_CONST) {
			if (isglobal) {
			    gvar[g] = sem.i;
			}
			else {
			    emit (OP_LOADC);
			    emit (sem.i);
			    emit (OP_I2F);
			    emit (OP_STORE);
			    emit (l+1);
			}
		    }
		    else {
			prerror ("incompatible types in initialization");
		    }
		}
		else {
		    prerror ("wrong type");
		}

		gettok ();		/* const */
	    }
	    else {			/* variable declaration type checking */
		if (type == TYPE_VOID) {
		    prerror ("storage size isn't known");
		}
	    }

	}
	if (tok == ';') {
	    gettok ();		/* ; */
	}      
    }
}


/* 
   expression
*/
static void expr (int etype, context_t *parent_cont) {
    context_t cont, cont1;
    int oldtype, optype;
    int tok1;
    int i, addr;
  
    switch (etype) {
    case 0:			/* expr */
	expr (1, &cont);
	*parent_cont = cont;
	break;
    case 1:			/* assign operator is right associative, 
				   so use recursion to implement right to left order */
	expr (2, &cont);
	if (tok == '=') {
	    if (!cont.islvalue) {
		prerror ("invalid lvalue");
	    }
	    cont1 = cont;
	    gettok (); /* = */
	    expr (1, &cont);
	    if ((cont1.type == TYPE_INT) && (cont.type == TYPE_FLOAT)) {
		emit (OP_F2I);
	    }
	    else if ((cont1.type == TYPE_FLOAT) && (cont.type == TYPE_INT)) {
		emit (OP_I2F);
	    }
	    if (cont1.isgvar) {
		emit (OP_STOREG);
		emit (cont1.gvar_addr);
		emit (OP_POP);
		emit (OP_LOADG);
		emit (cont1.gvar_addr);
	    }
	    else {
		emit (OP_STORE);
		emit (cont1.lvar_addr);
		emit (OP_POP);	/* pop old lvalue */
		emit (OP_LOAD);	/* push new lvalue */
		emit (cont1.lvar_addr);
	    }
	}
	*parent_cont = cont;
	break;
    case 2:			/* logical or */
	expr (3, &cont);
	*parent_cont = cont;
	while (tok == TOK_L_OR) {
	    parent_cont->islvalue = 0;
	    gettok ();
	    if (cont.type == TYPE_FLOAT) {
		emit (OP_F2I);
	    }
	    expr (3, &cont);
	    if (cont1.type == TYPE_FLOAT) {
		emit (OP_F2I);
	    }
	    emit (OP_L_OR);
	    parent_cont->type = TYPE_INT;
	}
	break;
    case 3:			/* logical and */
	expr (4, &cont);
	*parent_cont = cont;
	while (tok == TOK_L_AND) {
	    parent_cont->islvalue = 0;
	    gettok ();
	    if (cont.type = TYPE_FLOAT) {
		emit (OP_F2I);
	    }
	    expr (4, &cont);
	    if (cont.type = TYPE_FLOAT) {
		emit (OP_F2I);
	    }
	    emit (OP_L_AND);
	    parent_cont->type = TYPE_INT;
	}
	break;
    case 4:			/* equality */
	expr (5, &cont);
	oldtype = cont.type;
	*parent_cont = cont;
	while (tok == TOK_EQ || tok == TOK_NEQ) {
	    tok1 = tok;
	    parent_cont->islvalue = 0;
	    gettok ();
	    expr (5, &cont);
	    if (oldtype == TYPE_INT && cont.type == TYPE_FLOAT) {
		emit (OP_SWAP);
		emit (OP_I2F);
		optype = TYPE_FLOAT;
	    }
	    else if (oldtype == TYPE_FLOAT && cont.type == TYPE_INT) {
		emit (OP_I2F);
		optype = TYPE_FLOAT;
	    }
	    else if (oldtype == TYPE_FLOAT && cont.type == TYPE_FLOAT) {
		optype = TYPE_FLOAT;
	    }
	    else {
		optype = TYPE_INT;
	    }
	    oldtype = TYPE_INT;
	    if (optype == TYPE_FLOAT) {
		(tok == TOK_EQ) ? emit (OP_FEQ) : emit (OP_FNE);
	    }
	    else {
		(tok == TOK_EQ) ? emit (OP_IEQ) : emit (OP_INE);
	    }
	    parent_cont->type = TYPE_INT;
	}
	break;
    case 5:			/* relational */
	expr (6, &cont);
	*parent_cont = cont;
	oldtype = cont.type;
	while ((tok == '>') || (tok == '<') || (tok == TOK_GE) || (tok == TOK_LE)) {
	    parent_cont->islvalue = 0;
	    tok1 = tok;
	    gettok ();
	    expr (6, &cont);
	    if (oldtype == TYPE_INT && cont.type == TYPE_FLOAT) {
		emit (OP_SWAP);
		emit (OP_I2F);
		optype = TYPE_FLOAT;
	    }
	    else if (oldtype == TYPE_FLOAT && cont.type == TYPE_INT) {
		emit (OP_I2F);
		emit (OP_SWAP);
		optype = TYPE_FLOAT;
	    }
	    else if (oldtype == TYPE_FLOAT && cont.type == TYPE_FLOAT) {
		emit (OP_SWAP);
		optype = TYPE_FLOAT;
	    }
	    else {
		emit (OP_SWAP);
		optype = TYPE_INT;
	    }
	    oldtype = TYPE_INT;
	    switch (tok1) {
	    case '>':
		optype == TYPE_INT ? emit (OP_IGT) : emit (OP_FGT);
	    break;
	    case '<':
		optype == TYPE_INT ? emit (OP_ILT) : emit (OP_FLT);
	    break;
	    case TOK_GE:
		optype == TYPE_INT ? emit (OP_IGE) : emit (OP_FGE);
		break;
	    case TOK_LE:
		optype == TYPE_INT ? emit (OP_ILE) : emit (OP_FLE);
		break;
	    }
	    parent_cont->type = TYPE_INT;
	}
	break;
    case 6:			/* addictive */
	expr (7, &cont);
	*parent_cont = cont;
	oldtype = cont.type;
	while ((tok == '+') || (tok == '-')) {
	    if (tok == '+') {
		gettok ();
		expr (7, &cont);
		if (oldtype == TYPE_INT && cont.type == TYPE_FLOAT) {
		    emit (OP_SWAP);
		    emit (OP_I2F);
		    optype = TYPE_FLOAT;
		}
		else if (oldtype == TYPE_FLOAT && cont.type == TYPE_INT) {
		    emit (OP_I2F);
		    optype = TYPE_FLOAT;
		}
		else if (oldtype == TYPE_FLOAT && cont.type == TYPE_FLOAT) {
		    optype = TYPE_FLOAT;
		}
		else {
		    optype = TYPE_INT;
		}

		(optype == TYPE_FLOAT) ? emit (OP_FADD) : emit (OP_IADD);
	    }
	    else {
		gettok ();
		expr (7, &cont);
		if (oldtype == TYPE_INT && cont.type == TYPE_FLOAT) {
		    emit (OP_SWAP);
		    emit (OP_I2F);
		    optype = TYPE_FLOAT;
		}
		else if (oldtype == TYPE_FLOAT && cont.type == TYPE_INT) {
		    emit (OP_I2F);
		    emit (OP_SWAP);
		    optype = TYPE_FLOAT;
		}
		else if (oldtype == TYPE_FLOAT && cont.type == TYPE_FLOAT) {
		    emit (OP_SWAP);
		    optype = TYPE_FLOAT;
		}
		else {
		    emit (OP_SWAP);
		    optype = TYPE_INT;
		}
		optype == TYPE_FLOAT ? emit (OP_FSUB) : emit (OP_ISUB);
	    }
	    oldtype = optype;
	    parent_cont->type = oldtype;
	}
	break;
    case 7:			/* multiplicative */
	expr (8, &cont);
	*parent_cont = cont;
	oldtype = cont.type;
	while ((tok == '*') || (tok == '/') || (tok == '%')) {
	    parent_cont->islvalue = 0;
	    switch (tok) {
	    case '*':
		gettok ();
		expr (8, &cont);
		if (oldtype == TYPE_FLOAT && cont.type == TYPE_INT) {
		    emit (OP_I2F);
		    optype = TYPE_FLOAT;
		}
		else if (oldtype == TYPE_INT && cont.type == TYPE_FLOAT) {
		    emit (OP_SWAP);
		    emit (OP_I2F);
		    optype = TYPE_FLOAT;
		}
		else if (oldtype == TYPE_FLOAT && cont.type == TYPE_FLOAT) {
		    optype = TYPE_FLOAT;
		}
		else {
		    optype = TYPE_INT;
		}
		optype == TYPE_FLOAT ? emit (OP_FMUL) : emit (OP_IMUL);
		oldtype = optype;
		break;
	    case '/':
		gettok ();
		expr (8, &cont);
		if (oldtype == TYPE_INT && cont.type == TYPE_FLOAT) {
		    emit (OP_SWAP);
		    emit (OP_I2F);
		    optype = TYPE_FLOAT;
		}
		else if (oldtype == TYPE_FLOAT && cont.type == TYPE_INT) {
		    emit (OP_I2F);
		    emit (OP_SWAP);
		    optype = TYPE_FLOAT;
		}
		else if (oldtype == TYPE_FLOAT && cont.type == TYPE_FLOAT) {
		    emit (OP_SWAP);
		    optype = TYPE_FLOAT;
		}
		else {
		    emit (OP_SWAP);
		    optype = TYPE_INT;
		}
		optype == TYPE_FLOAT ? emit (OP_FDIV) : emit (OP_IDIV);
		oldtype = optype;
		break;
	    case '%':
		gettok ();
		expr (8, &cont);
		if (oldtype == TYPE_FLOAT || cont.type == TYPE_FLOAT) {
		    prerror ("invalid operands to binary %");
		}
		emit (OP_SWAP);
		emit (OP_MOD);
		oldtype = TYPE_INT;
		optype = TYPE_INT;
		break;
	    }
	    oldtype = optype;
	    parent_cont->type = optype;
	}
	break;
    case 8:			/* unary */
	if (tok == '+') {
	    parent_cont->islvalue = 0;
	    gettok ();
	    expr (9, &cont);
	}
	else if (tok == '-') {
	    parent_cont->islvalue = 0;
	    gettok ();
	    expr (9, &cont);
	    emit (OP_NEG);
	}
	else if (tok == '!') {
	    parent_cont->islvalue = 0;
	    gettok ();
	    expr (9, &cont);
	    emit (OP_NOT);
	}
	else {
	    expr (9, &cont);
	    *parent_cont = cont;
	}
	break;
    case 9:			/* primary */
	if (tok == TOK_ID) {
	    if ((addr = findlsym (sem.id)) != -1) {
		parent_cont->isgvar = 0;
		parent_cont->lvar_addr = lvarcont[addr].lvar_addr;
	    }
	    else if ((addr = findgsym (sem.id)) != -1) {
		parent_cont->isgvar = 1;
		parent_cont->gvar_addr = addr;
	    }
	    else {
		prerror ("undeclared variable");
	    }
	    parent_cont->islvalue = 1;
	    gettok ();
	    if (tok == '(') {
		parent_cont->islvalue = 0;
		gettok ();
		i = 0;
		if (parent_cont->isgvar) {
		    if (gvarcont[addr].isfunc) {
			while (tok != ')') {
			    if (tok == ',') {
				gettok ();
			    }
			    oldtype = gvarcont[addr].argtype[i];
			    expr (0, &cont);
			    if (oldtype == TYPE_INT && cont.type == TYPE_FLOAT) {
				emit (OP_F2I);
			    }
			    else if (oldtype == TYPE_FLOAT && cont.type == TYPE_INT) {
				emit (OP_I2F);
			    }
			    i++;
			}
			if (i > gvarcont[addr].numargs) {
			    prerror ("too many arguments");
			}
			if (i < gvarcont[addr].numargs) {
			    prerror ("too few arguments");
			}
		    }
		    else {
			prerror ("not a function");
		    }
		}
		else {
		    prerror ("no local function");
		}
		gettok ();		/* ) */
		if (gvarcont[addr].isbuildin) {
		    emit (OP_CALLN);
		    emit (addr);
		}
		else {
		    emit (OP_CALL);
		    if (gvarcont[addr].isdecl) {
			emit (0);
			funcpatch (addr, numops-1);
		    }
		    else {
			emit (gvar[addr]);
		    }
		}
		emit (OP_COMPRESS);
		emit (gvarcont[addr].numargs);	/* return hack */

		parent_cont->type=gvarcont[addr].type;
	    }
	    else {
		if (parent_cont->isgvar) {
		    parent_cont->type = gvarcont[addr].type;
		    emit (OP_LOADG);
		    emit (parent_cont->gvar_addr);
		}
		else {
		    parent_cont->type = lvarcont[addr].type;
		    emit (OP_LOAD);
		    emit (parent_cont->lvar_addr);
		}
	    }


	}
	else if (tok == '(') {
	    gettok ();
	    expr (0, &cont);
	    *parent_cont = cont;
	    expect (')');
	}
	else {
	    expr (10, &cont);
	    parent_cont->type = cont.type;
	    parent_cont->islvalue = 0;
	}
	break;

    case 10:			/* const */
	parent_cont->islvalue = 0;
	switch (tok) {
	case TOK_INT_CONST:
	    parent_cont->type = TYPE_INT;
	    emit (OP_LOADC);
	    emit (sem.i);
	    break;
	case TOK_FLOAT_CONST:
	    parent_cont->type = TYPE_FLOAT;
	    emit (OP_LOADC);
	    emit (sem.i);
	    break;
	case TOK_STRING_CONST:
	    parent_cont->type = TYPE_INT;
	    emit (OP_LOADC);
	    emit (addstr (sem.str));
	    break;
	default:
	    prerror ("unknown const type");
	}
	gettok ();
	break;
    }
}

enum {JMP_BREAK, JMP_CONTINUE};
typedef struct {
    int opaddr;
    int type;
} jmp_t;
jmp_t jmptab[10000];
int numjmps;
int jmpstart;
int jmpstack[100];
int jmptos;			/* top of jmp stack */
/* 
 */
static void beginjmp () {
    jmpstart += numjmps;
    jmpstack[jmptos++] = numjmps;
    numjmps = 0;
}

static void endjmp (int b, int c) {
    jmp_t *start;
    int i;
  
    start = jmptab + jmpstart;
    for (i=0; i<numjmps; i++) {
	switch (start[i].type) {
	case JMP_BREAK:
	    bytecode[start[i].opaddr] = b;
	    break;
	case JMP_CONTINUE:
	    bytecode[start[i].opaddr] = c;
	    break;
	}
    }
    numjmps = jmpstack[--jmptos];
    jmpstart -= numjmps;
}

static void addjmp (int opaddr, int type) {
    jmptab[jmpstart+numjmps].opaddr = opaddr;
    jmptab[jmpstart+numjmps++].type = type;
}


/* 
   what is block? everything is block!
*/
static void block (int isfuncdecl, int funcaddr, int breakable) {
    context_t cont;
    int lab, lab1, lab2, lab3;
    int brkaddr, cntaddr;
    int i;

    if (tok == '{') {
	gettok ();
	if (isfuncdecl) {
	    emit (OP_ENLARGE);
	    emit (0);
	    lab = genlab ();
	    patchlab (lab, numops-1);
	    decl (0);
	    gvarcont[funcaddr].numlvars = numlsyms - gvarcont[funcaddr].numargs;
	    backpatchlab (lab, gvarcont[funcaddr].numlvars);
	    for (i=0; i<gvarcont[funcaddr].numlvars; i++) { /* addressing part 2: local vars */
		lvarcont[i + gvarcont[funcaddr].numargs].lvar_addr = i+1;
	    }
	    while (tok != '}') {
		block (0, funcaddr, breakable);
	    }
	    if (gvarcont[funcaddr].type == TYPE_VOID) {
		emit (OP_RET);
	    }
	    else {
		emit (OP_RETN);
	    }            
	    gettok ();  
/*      resetlsym ();          */
	}
	else {
	    while (tok != '}') {
		block (0, funcaddr, breakable);
	    }
	    gettok ();
	}      
    }
    else if (tok == TOK_IF) {
	gettok ();
	expect ('(');
	expr (0, &cont);
	expect (')');
	lab = genlab ();
	emit (OP_JZ);		/* if false, jump to the end addr of if statement */
	emit (0);
	patchlab (lab, numops-1);	/* but we don't know the addr yet */
	block (0, funcaddr, breakable);
	if (tok == TOK_ELSE) {
	    gettok ();
	    emit (OP_JMP);		/* if there's else statement, then the if body should jump over that */
	    emit (0);
	    backpatchlab (lab, numops);	/* the addr of if end in this situation */
	    lab1 = genlab ();
	    patchlab (lab1, numops-1);
	    block (0, funcaddr, breakable);
	    backpatchlab (lab1, numops);
	}
	else {
	    backpatchlab (lab, numops);	/* now we know the addr of the end of if statement*/
	}
    }
    else if (tok == TOK_WHILE) {
	gettok ();
	expect ('(');
	beginjmp ();    
	cntaddr = numops;    
	lab1 = numops;
	expr (0, &cont);
	expect (')');
	lab = genlab ();
	emit (OP_JZ);
	emit (0);
	patchlab (lab, numops-1);
	block (0, funcaddr, 1);
	emit (OP_JMP);
	emit (lab1);
	backpatchlab (lab, numops);
	brkaddr = numops;
	endjmp (brkaddr, cntaddr);
    }
    else if (tok == TOK_FOR) {
	gettok ();
	expect ('(');
	if (tok != ';') {
	    expr (0, &cont);		/* first expr */
	}
	expect (';');
	if (tok != ';') {
	    lab3 = numops;
	    expr (0, &cont);		/* second expr */
	    lab = genlab ();
	    emit (OP_JZ);
	    emit (0);
	    patchlab (lab, numops-1);/* if expr == false, jump to the end of while statement, but now we don't know the addr */
	    emit (OP_JMP);
	    emit (0);
	    lab2 = genlab ();
	    patchlab (lab2, numops-1); /* jump over the third expr */
	}
	expect (';');
	if (tok != ')') {
	    lab1 = numops;
	    cntaddr = numops;
	    expr (0, &cont);		/* third expr */
	    emit (OP_JMP);		/* jump to the 2nd expr */
	    emit (lab3);
	    backpatchlab (lab2, numops);	
	}
	expect (')');
	block (0, funcaddr, 1);
	emit (OP_JMP);
	emit (lab1);
	brkaddr = numops;
	backpatchlab (lab, numops); /* now we know the addr of the end of while statement */
	endjmp (brkaddr, cntaddr);
    }
    else if (tok == TOK_BREAK) {
	if (!breakable) {
	    prerror ("break statement not within loop");
	}
	gettok ();
	emit (OP_JMP);
	emit (0);
	addjmp (numops-1, JMP_BREAK);
    }
    else if (tok == TOK_CONTINUE) {
	if (!breakable) {
	    prerror ("continue statement not within loop");
	}
	gettok ();
	emit (OP_JMP);
	emit (0);
	addjmp (numops-1, JMP_CONTINUE);
    }
    else if (tok == TOK_RETURN) {
	gettok ();
	if (tok != ';') {
	    expr (0, &cont);
	    emit (OP_RETN);
	}
	else {
	    emit (OP_RET);
	}
	expect (';');
    }
    else {			/* expression */
	if (tok != ';') {
	    expr (0, &cont);
	    emit (OP_POP);
	}
	expect (';');
    }
}



/* 
 */
void compilefile(char *filename) {
    int i;
  
    f = fopen (filename, "r");
    if (!f) {
	exit (1);
    }    
    i = setjmp (jb);

    if (!i) {
	getch ();
	gettok ();
	emit (OP_END);		/* quick hack: exit point is the first word of the bytecode  */
	decl (1);
	funcbackpatch ();
    }
    else {
	fclose (f);
	exit (1);    
    }
  
}




void beginparse(char *filename)
{
    int i;

    f = fopen(filename, "r");
    i = setjmp(jb);

    getch();
}


void endparse()
{
    fclose(f);
}

/* 
 */
void savebytecode() {
}


/* 
 */
void loadbytecode() {
}

int readhash(char *key) {
}

void writehash(char *key, int value) {
}


/* 
 */
static int findfunc (char *funcname) {
    int i;

    for (i=0; i<numgsyms; i++) {
	if (!gvarcont[i].isfunc) {
	    continue;
	}
	if (!strcmp (gsym[i], funcname)) {
	    return gvar[i];
	}
    }
    return -1;
}

/* 
 */
int execbytecode(char *funcname, int numargs, int arg[]) {
    int tmp;
    int op;
    float f, f1;
    int i;

    pc = findfunc (funcname);
    if (pc == -1) {
	printf ("can't find function \"%s\"", funcname);
	return;
    }
    for (i=0; i<numargs; i++) {
	mem[i] = arg[i];
    }
    mem[numargs] = 0;
    mem[numargs+1] = 0;
    sp = numargs+1;
    fp = sp;
    while ((op = bytecode[pc]) != OP_END) {
	switch (op) {
	case OP_POP:
	    sp --;
	    pc ++;
	    break;
	case OP_ENLARGE:
	    op = bytecode[++pc];
	    sp += op;
	    pc ++;
	    break;
	case OP_COMPRESS:
	    op = bytecode[++pc];
	    mem[sp-op] = mem[sp];
	    sp -= op;
	    pc ++;
	    break;
	case OP_SWAP:
	    tmp = mem[sp];
	    mem[sp] = mem[sp-1];
	    mem[sp-1] = tmp;
	    pc ++;
	    break;
	case OP_DUP:
	    mem[sp+1] = mem[sp];
	    sp ++;
	    pc ++;
	    break;
	case OP_LOADC: 
	    op = bytecode[++pc]; 
	    mem[++sp] = op; 
	    pc ++;
	    break;
	case OP_LOAD:
	    op = bytecode[++pc];
	    mem[++sp] = mem[fp + op];
	    pc++;
	    break;
	case OP_STORE:
	    op = bytecode[++pc];
	    mem[fp + op] = mem[sp--];
	    pc ++;
	    break;
	case OP_LOADG:
	    op = bytecode[++pc];
	    mem[++sp] = gvar[op];
	    pc ++;
	    break;
	case OP_STOREG:
	    op = bytecode[++pc];
	    gvar[op] = mem[sp--];
	    pc ++;
	    break;
	case OP_IADD:
	    mem[sp-1] = mem[sp] + mem[sp-1];
	    sp --;
	    pc ++;
	    break;
	case OP_ISUB:
	    mem[sp-1] = mem[sp] - mem[sp-1];
	    sp --;
	    pc ++;
	    break;
	case OP_IMUL:
	    mem[sp-1] = mem[sp] * mem[sp-1];
	    sp --;
	    pc ++;
	    break;
	case OP_IDIV:
	    mem[sp-1] = mem[sp] / mem[sp-1];
	    sp --;
	    pc ++;
	    break;
	case OP_MOD:
	    mem[sp-1] = mem[sp] % mem[sp-1];
	    sp --;
	    pc ++;
	    break;
	case OP_L_OR:
	    mem[sp-1] = mem[sp] || mem[sp-1];
	    sp --;
	    pc ++;
	    break;
	case OP_L_AND:
	    mem[sp-1] = mem[sp] && mem[sp-1];
	    sp --;
	    pc ++;
	case OP_IGT:			/* greater than */
	    mem[sp-1] = mem[sp] > mem[sp-1];
	    sp --;
	    pc ++;
	    break;
	case OP_ILT:			/* less than */
	    mem[sp-1] = mem[sp] < mem[sp-1];
	    sp --;
	    pc ++;
	    break;
	case OP_IGE:			/* greater equal */
	    mem[sp-1] = mem[sp] >= mem[sp-1];
	    sp --;
	    pc ++;
	    break;
	case OP_ILE:		/* less equal */
	    mem[sp-1] = mem[sp] <= mem[sp-1];
	    sp --;
	    pc ++;
	    break;
	case OP_IEQ:		/* equal */
	    mem[sp-1] = mem[sp] == mem[sp-1];
	    sp --;
	    pc ++;
	    break;
	case OP_INE:		/* not equal */
	    mem[sp-1] = mem[sp] != mem[sp-1];
	    sp --;
	    pc ++;
	    break;
	case OP_I2F:
	    f = mem[sp];
	    mem[sp] = *((int *) &f);
	    pc ++;
	    break;
	case OP_F2I:
	    f = *((float *) &(mem[sp]));
	    mem[sp] = f;
	    pc ++;
	    break;
	case OP_FADD:
	    f = *((float *) &mem[sp]);
	    f1 = *((float *) &mem[--sp]);
	    f1 = f + f1;
	    mem[sp] = *((int *) &f1);
	    pc ++;
	    break;
	case OP_FSUB:
	    f = *((float *) &mem[sp]);
	    f1 = *((float *) &mem[--sp]);
	    f1 = f - f1;
	    mem[sp] = *((int *) &f1);
	    pc ++;
	    break;
	case OP_FMUL:
	    f = *((float *) &mem[sp]);
	    f1 = *((float *) &mem[--sp]);
	    f1 = f * f1;
	    mem[sp] = *((int *) &f1);
	    pc ++;
	    break;
	case OP_FDIV:
	    f = *((float *) &mem[sp]);
	    f1 = *((float *) &mem[--sp]);
	    f1 = f / f1;
	    mem[sp] = *((int *) &f1);
	    pc ++;
	    break;
	case OP_FGT:			/* greater than */
	    f = *((float *) &mem[sp]);
	    f1 = *((float *) &mem[--sp]);
	    f1 = f > f1;
	    mem[sp] = *((int *) &f1);
	    pc ++;
	    break;
	case OP_FLT:			/* less than */
	    f = *((float *) &mem[sp]);
	    f1 = *((float *) &mem[--sp]);
	    f1 = f < f1;
	    mem[sp] = *((int *) &f1);
	    pc ++;
	    break;
	case OP_FGE:			/* greater equal */
	    f = *((float *) &mem[sp]);
	    f1 = *((float *) &mem[--sp]);
	    f1 = f >= f1;
	    mem[sp] = *((int *) &f1);
	    pc ++;
	    break;
	case OP_FLE:		/* less equal */
	    f = *((float *) &mem[sp]);
	    f1 = *((float *) &mem[--sp]);
	    f1 = f <= f1;
	    mem[sp] = *((int *) &f1);
	    pc ++;
	    break;
	case OP_FEQ:		/* equal */
	    f = *((float *) &mem[sp]);
	    f1 = *((float *) &mem[--sp]);
	    f1 = f == f1;
	    mem[sp] = *((int *) &f1);
	    pc ++;
	    break;
	case OP_FNE:		/* not equal */
	    f = *((float *) &mem[sp]);
	    f1 = *((float *) &mem[--sp]);
	    f1 = f != f1;
	    mem[sp] = *((int *) &f1);
	    pc ++;
	    break;
	case OP_JMP:
	    pc = bytecode[pc + 1];
	    break;
	case OP_JZ:
	    op = bytecode[++pc];
	    if (!mem[sp--]) {
		pc = op;
	    }
	    else {
		pc ++;
	    }
	    break;
	case OP_JNZ:
	    op = bytecode[++pc];
	    if (mem[sp--]) {		/* pop */
		pc = op;
	    }
	    else {
		pc ++;
	    }
	    break;
	case OP_CALL:
	    op = bytecode[++pc];
	    mem[++sp] = pc + 1;
	    pc = op;
	    mem[++sp] = fp;
	    fp = sp;
	    break;
	case OP_CALLN:		/* call native functions */
	    op = bytecode[++pc];
	    if (gvarcont[op].type == TYPE_VOID) {
		mem[sp+1] = ((buildinfp_t ) gvar[op]) (sp);
		sp ++;
	    }
	    else {
		((buildinfpvoid_t) gvar[op]) (sp);
		mem[++sp] = 0;
	    }
	    pc ++;
	    break;
	case OP_RET:
	    sp = fp-1;
	    pc = mem[fp - 1];
	    fp = mem[fp];
	    mem[sp] = 0;
	    break;
	case OP_RETN:
	    pc = mem[fp-1];
	    mem[fp-1] = mem[sp];
	    sp = fp-1;
	    fp = mem[fp];
	    break;
	case OP_NOP:
	    pc ++;
	    break;
	}
    }
    return mem[sp];
}



#if 0

int main(int argc, char *argv[])
{
    if(argc==1){
	printf("usage: lc <filename>\n");
	return 0;
    }


    printf("*****************************start compiling******************************\n\n");
    compilefile(argv[1]);
    printf("*********************************complete*********************************\n\n");


    execbytecode("main",0,NULL);

    return 0;
}

#endif
