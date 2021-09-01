/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include"smoothness.h"
#define FP_SIZE 2
 
static int calculate_newton_indices(int *ind, int k)
{
    if (k == 1)
        return 0;
    ind[0] = k;
    return 1+calculate_newton_indices(ind+1, (1+k)/2);
}


/*
 Iteration formula:
    g <- ( 2*g*t^(deg(g)+deg(f)) - f*g^2 ) div t^xxx
*/


void GF2X_msb_preinverse(GF2X& inv_f, GF2X& f, int k)
{
    ASSERT(k < 1024);

    int indices[10]; // 2^10 = 1024. More than enough.
    int l = calculate_newton_indices(indices, k+1);
    ASSERT(indices[l-1] == 2);
    ASSERT(indices[0] == k+1);
    GF2X ff, g, tmp;

    GF2X_set_ti(g, 0);
    for (int i = l-1; i >= 0; --i) 
    {
        int ind = indices[i];
       
        sqr(tmp, g);
        
        RightShift(ff, f, max(0, deg(f)-ind)); 
        mul(g, tmp, ff);    
        RightShift(g, g, deg(g)-ind);
    }
    RightShift(inv_f,g, deg(g)-k);
    ff.kill();
    g.kill();
    tmp.kill();
}

void GF2X_rem_precomp(GF2X& r, GF2X& pq, 
        GF2X& m, GF2X& invm)
{
    GF2X pqt;

    int d = deg(m);
    ASSERT_ALWAYS(deg(pq) <= deg(m) + deg(invm));
    //GF2X_init(&pqt);
    RightShift(pqt, pq, d);
    mul(pqt, pqt, invm);  // should be a mul-high
    RightShift(pqt,pqt, deg(invm));
    mul(pqt,pqt, m);     // should be a mul-low
    sub(r, pq, pqt);
    ASSERT_ALWAYS(deg(r) < deg(m));
    pqt.kill();
}


// The return value is 
//   - 0 if P is not smooth
//   - an upper bound on the degree of the largest factor otherwise
// NB: this algorithm is allowed to fail. If there is a factor of even
// multiplicity that has a degree > B, it is not detected.

int GF2X_is_smooth(GF2X PP, int B)
{
    int B2 = (1+B)>>1;  // Ceiling(B/2)

    if (deg(PP) <= B)
        return deg(PP);

    GF2X P;  // monic version of PP
    {
        GF2 lc;
       lc=coeff(PP, deg(PP));
        if (IsOne(lc))
            P=PP;
        else
            div(P, PP, lc);
    }
    GF2X dP;
    diff(dP, P);
    
    GF2X tqi, acc, t, tmp;
    GF2X invP2, invPq;

    GF2X_msb_preinverse(invP2, P, deg(P)-2);
    invPq=invP2;

    GF2X_set_ti(t, 1);

    int q = FP_SIZE;
    int qi = q;
    int i = 1;
    while (qi < deg(P)) {
        qi *= q;
        i++;
    }
    GF2X_set_ti(tqi, qi);
    GF2X_rem_precomp(tqi, tqi, P, invPq);

    while (i < B2) {
        sqr(tqi, tqi);
        GF2X_rem_precomp(tqi, tqi, P, invPq);
        i++;
    }

    int smooth = 0;
    sub(acc, tqi, t);
    mul(acc, acc, dP);
    GF2X_rem_precomp(acc, acc, P, invP2);
    while (i < B) {
        sqr(tqi, tqi);
        GF2X_rem_precomp(tqi,tqi,P,invPq);
        i++;
        sub(tmp,tqi,t);
        mul(acc, acc, tmp);
        GF2X_rem_precomp(acc,acc,P,invP2);
        if (IsZero(acc)) {
            smooth = 1;
            break;
        }
    }

    P.kill();
    dP.kill();
    tqi.kill();
    t.kill();
    tmp.kill();
    acc.kill();
    invP2.kill();
    invPq.kill();

    if (!smooth)
        return 0;
    else
        return i;
}
