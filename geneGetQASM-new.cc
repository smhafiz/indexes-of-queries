#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

using namespace std;

#define LO 			0
#define HI 			1

#define CARRY 			1
#define MUL			2
#define MADC 			3

#define CARRY_IN(i)		(state[i] |= CARRY)
#define MUL_IN(i)		(state[i] |= MUL)

#define IS_OCCUPIED(i)		(state[i] & MADC ? true : false)
#define DO_CARRY_OUT(i)		(state[i] & MUL ? true : false)
#define CARRY_OUT_FLAG(i)	(DO_CARRY_OUT(i) ? ".cc" : "")
#define CARRY_IN_FLAG(b)	(b ? "c" : "")
#define LO_OR_HI(p)		(p.first == HI ? ".hi" : ".lo")

#define GET_DEST_REG(i)  	(i < n-1 ? i+1 : i-n+1)//(i < n-1 ? n+i-1 : i-1)
#define GET_FIRST_REG(p)	(n+1+p.second)//(2*n-1+p.second)
#define GET_SECOND_REG(i,p)	(2*n+1+i-p.second-p.first)//(3*n-1+i-p.second-p.first)

vector< vector< pair<int,int> > > *pairs;
int * state;
char * r0 = new char[8];
char * r1 = new char[8];
char * r2 = new char[8];

void propagate(int i, int n, bool c)
{
    bool cc = DO_CARRY_OUT(i);
    sprintf(r0, "%%%u", GET_DEST_REG(i));
    //if (!c) cout << "\n";
    if ((*pairs)[i].empty())
    {
	printf("\t\"add%s%s.u32\t%3s,%3s,  0    ;\"\\n\\t", CARRY_IN_FLAG(c), CARRY_OUT_FLAG(i), r0, (IS_OCCUPIED(i) ? r0 : "  0"));
	sprintf(r0, "r%u", GET_DEST_REG(i));
	printf("\t//%3s%s= c", r0, (IS_OCCUPIED(i) ? "+" : " ")); 
	CARRY_IN(i);
    }
    else
    {
	auto p = (*pairs)[i].back();
	(*pairs)[i].pop_back();
	sprintf(r1, "%%%u", GET_FIRST_REG(p));
	sprintf(r2, "%%%u", GET_SECOND_REG(i,p));
	printf("\t\"mad%s%s%s.u32\t%3s,%3s,%3s,%3s;\"\\n\\t", CARRY_IN_FLAG(c), LO_OR_HI(p), CARRY_OUT_FLAG(i), r0, r1, r2, (IS_OCCUPIED(i) ? r0 : "  0"));
	sprintf(r0, "r%u", GET_DEST_REG(i));
	sprintf(r1, "r%u", GET_FIRST_REG(p));
	sprintf(r2, "r%u", GET_SECOND_REG(i,p));
	printf("\t//%3s%s=[%3s*%3s]%s%s", r0, (IS_OCCUPIED(i) ? "+" : " "), r1, r2, LO_OR_HI(p), (c ? "+c" : "  "));
	MUL_IN(i);
    }
    if (i < n - 1) cout << "  (r" << (i - n + 1) << " => r" << GET_DEST_REG(n + i) << ")";
    cout << "\n";
    if (cc && i < 2 * n - 1) propagate(i + 1, n, true);
}

int main(int argc, char ** argv)
{
    int n = atoi(argv[1]);
    state = new int[2 * n];
    pairs = new vector< vector< pair<int,int> > >();
    pairs->resize(2 * n);

    for (int i = 0; i < 2 * n; ++i)
    {
	for (int j = 0; j <= i; j++)
	{
	    if (i - j < n - 1 && j < n)
	    {
		if (i == 0 && j == 0) continue;
		(*pairs)[i].push_back(make_pair(LO, j));
		(*pairs)[i+1].push_back(make_pair(HI, j));
	    }
	}
    }

    cout << "__device__ __forceinline__ uint" << (32*n) << " get_q(const uint" << (32*n) << " & b, const uint" << (32*(n-1)) << " & c)\n";
    cout << "{\n";
    cout << "    uint tmp;\n";
    cout << "    uint" << (32*n) << " a;\n";
    cout << "    uint * _a = (uint *) &a;\n";
    cout << "    uint * _b = (uint *) &b;\n";
    cout << "    uint * _c = (uint *) &c;\n";
    auto p = make_pair(1,0);
    sprintf(r0, "%%%u", GET_DEST_REG(n+1));
    sprintf(r1, "%%%u", GET_FIRST_REG(p));
    sprintf(r2, "%%%u", GET_SECOND_REG(1,p));
    printf("    asm(\"mul.hi.u32\t%3s,%3s,%3s    ;\"\\n\\t", r0, r1, r2);
    sprintf(r0, "r%u", GET_DEST_REG(n+1));
    sprintf(r1, "r%u", GET_FIRST_REG(p));
    sprintf(r2, "r%u", GET_SECOND_REG(1,p));
    printf("\t//%3s =[%3s*%3s].hi  ", r0, r1, r2);
    cout << "  (r" << (2-n) << " => r" << GET_DEST_REG(n+1) <<  ")\n";
    MUL_IN(1);

    for (int i = 0; i < 2 * n; i++)
    {
	while (!(*pairs)[i].empty())
	{
	    propagate(i, n, false);
	}
    }

    // handle the implicit 1-word
    sprintf(r0, "%%%u", GET_DEST_REG(n));
    sprintf(r1, "%%%u", 2 * n - 1);
    printf("\t\"add.cc.u32\t%3s,%3s,%3s    ;\"\\n\\t",r0,r0,r1);
    sprintf(r0, "r%u", GET_DEST_REG(n));
    sprintf(r1, "r%u", 2 * n - 1);
    printf("\t//%3s+=%3s\n", r0, r1);
    for (int i = n; i < 2 * n - 2; i++)
    {
	sprintf(r0, "%%%u", GET_DEST_REG(i+1));
	sprintf(r1, "%%%u", i + n);
	printf("\t\"addc.cc.u32\t%3s,%3s,%3s    ;\"\\n\\t",r0,r0,r1);
	sprintf(r0, "r%u", GET_DEST_REG(i+1));
	sprintf(r1, "r%u", i + n);
	printf("\t//%3s+=%3s+c\n", r0, r1);
    }
    sprintf(r0, "%%%u", GET_DEST_REG(2 * n - 1));
    sprintf(r1, "%%%u", 3 * n - 2);
    printf("\t\"addc.u32\t%3s,%3s,%3s    ;\"\\n\\t",r0,r0,r1);
    sprintf(r0, "r%u", GET_DEST_REG(2 * n - 1));
    sprintf(r1, "r%u", 3 * n - 2);
    printf("\t//%3s+=%3s+c\n", r0, r1);

    sprintf(r0, "%%%u", GET_DEST_REG(2 * n - 1));
    sprintf(r1, "%%%u", 3 * n - 2);
    printf("  //\t\"addc.cc.u32\t%3s,%3s,%3s    ;\"\\n\\t\t// propagate overflow?\n",r0,r0,r1);
    sprintf(r0, "r%u", GET_DEST_REG(2 * n - 1));
    sprintf(r1, "r%u", 3 * n - 2);
    printf("  //\t\"addc.u32\t %%0,  0,  0    ;\"\\n\\t\t// overflow?\n");

    cout << "\t: \"=r\"(tmp)";
    for (int i = 0; i < n; i++) printf(", \"=r\"(_a[%u])", i);
    cout << "\n\t: \"r\"(_b[0])";
    for (int i = 1; i < n; i++) printf(", \"r\"(_b[%u])", i);
    for (int i = 0; i < n-1; i++) printf(", \"r\"(_c[%u])", i);
    cout << ");\n\n";
    cout << "    return a;\n";
    cout << "}\n\n";

    return 0;
}
