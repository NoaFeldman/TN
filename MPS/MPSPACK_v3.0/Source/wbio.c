#ifndef __WB_WBIO_C__
#define __WB_WBIO_C__

// ----------------------------------------------------------------- //
// ----------------------------------------------------------------- //

void WbPrintMatrixC (
    const wbMatrix<wbcomplex> &M, const char *istr, int space
){
    unsigned i,j; char dstr[64], sfmt[64];
    
    int rflag, iflag;

    printf ("%s%s", istr, istr[0] ? "\n" : "");

    if (space<1) sprintf (sfmt, " %%s");
    else         sprintf (sfmt, " %%%ds", space);

    for (i=0; i<M.dim1; i++) {
        for (j=0; j<M.dim2; j++) {
            
            rflag = std::isnan(M(i,j).r);
            if (!rflag) {
                if (M(i,j).r!=0.) {
                    rflag = (fabs(M(i,j).r)==1.) ? 1 : 2;
                    if (M(i,j).r<0) rflag = -rflag;
                }
            }

            iflag = std::isnan(M(i,j).i);
            if (!iflag) {
                if (M(i,j).i!=0.) {
                    iflag = (fabs(M(i,j).i)==1.) ? 1 : 2;
                    if (M(i,j).i<0) iflag = -iflag;
                }
            }


            if (rflag) {
                sprintf (dstr, "% g", M(i,j).r);
                if (iflag) {
                    if (abs(iflag)>1)
                         sprintf (dstr, "%s%+gi", dstr, M(i,j).i);
                    else sprintf (dstr, "%s%ci",  dstr, iflag==1 ? '+' : '-');
                }
            }
            else {
                if (iflag)
                    if (abs(iflag)>1)
                         sprintf (dstr, "% gi", M(i,j).i);
                    else sprintf (dstr, "%ci",  iflag==1 ? ' ' : '-');
                else
                    sprintf (dstr, "%g", M(i,j).r);
            }
            printf (sfmt, dstr);
        }
        printf ("\n");
    }
    printf ("\n");
}


wbvector<unsigned> Str2Idx(
    const char* s, unsigned offset, bool exitflag) {

    wbvector<unsigned> idx;

    int e=Str2Idx(s,idx,offset);
    if (e<=0) { e=-e;
       if (!s || unsigned(e)>=strlen(s))
          wblog(FL, "ERR %F: invalid index '%s'%N", s?s:"(null)");
       else {
          wbstring ss(' ', strlen(s)); ss[e]='^';
          wblog(FL, "ERR %F: invalid index '%s' (%d)%N%60s%s%N",
          s, e, "", ss.data);
       }
       if (exitflag)
       throw (char*)"err::invalid string index";
    }

    return idx;
};


template<class T>
int Str2Idx(const char* s, wbvector<T> &idx, unsigned offset) {

    unsigned i=0, k=0, l=0; int e=0;
    char c, bflag=0, sflag=1;

    if (s) {
       while ((c=s[l++])) {
          if (c>='0' && c<='9') { if (sflag) { ++k; sflag=0; }}
          else { sflag=1; }
       }
    }

    idx.init(k); k=-1; sflag=1;

    if (l<2) { return l; }; --l;

    T *d=idx.data;


    for (i=0; i<l; ++i) { c=s[i];
       if (c==' ') { sflag|=1; continue; }
       if (c==',') {
          if (sflag & 2) { e=i+1; break; }
          else { sflag|=2; continue; }
       }
       if (c>='0' && c<='9') {
          if (sflag) {
             d[++k]=c-'0'; sflag=0;
             if (k && d[k-1]<offset) { e=i+1; break; }
          }
          else d[k] = 10*d[k] + (c-'0');
          continue;
       }
       if (c=='[') {
          if ((++bflag)!=1 || int(k)>=0 || sflag!=1) { e=i+1; break; }
       } else
       if (c==']') {
          if ((++bflag)!=2 || sflag>1) e=i+1;
          for (++i; i<l; ++i) if (s[i]!=' ') { e=i+1; break; }
          if (e) break;
       } else { e=i+1; break; }
    }

    if (!sflag) { if (d[k++]<offset) e=i+1; } else
    if (sflag &2) e=i+1; else
    if (bflag==1) e=i+1; else
    if (k!=idx.len || i<l) e=i+1;

    if (offset) { for (i=0; i<k; ++i) d[i]-=offset; }

    if (e) return -e;

    return l;
};


#endif

