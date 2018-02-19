#ifndef __WB_LOG_C__
#define __WB_LOG_C__

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

void usage(const char *F, int L, const char* estr) {

    printf("\n%s\n", USAGE);
    if (estr[0]) wberror(F,L,estr);
}


void banner(unsigned n, const char *s, const char *istr, const char *fstr){

   unsigned i,r=0,l=0,m=0, s1=0, s2=0;

   if (s && s[0]) r=strlen(s);
   else wblog(FL,"ERR cannot repeat empty string!");

   if (istr && istr[0]) { s1=strlen(istr); printf("%s",istr); }
   if (fstr && fstr[0]) { s2=strlen(fstr); }

   m=s1+s2; if (n>=m) l=(n-m)/r;
   else wblog(FL,"ERR overall length too short (%d/%d)",n,m);

   for (i=0; i<l; i++) printf("%s",s);
   for (i=l*r+s1, l=n-s2; i<l; i++) printf("%s",s);

   if (fstr && fstr[0]) { printf("%s",fstr); }
   printf("\n");
};


void wberror(const char *F, int L, const char* istr) {

   printf("\n%s ERR\n\n%s\n\n", shortFL(F,L), istr);
   ExitMsg();
}

void wberror(
   const char *F, int L,
   const char *istr1,
   const char *istr2
){
   wblog(F,L,"ERR %s%N%N%s%N", istr1, istr2);
   ExitMsg();
}


inline void init_header(
   char *hstr, const char *file, int line,
   const char *time_stamp, const char *tag, int lenFL=21
){
   int i,k,n;

   if (file) {
      for (k=i=0; file[i]; i++) { if (file[i]=='/') k=i+1; }

      snprintf(hstr,120,"%.100s:%d", file+k, line);

      n=strlen(hstr); if (n>lenFL) {
          i=lenFL/2; hstr[i]=hstr[i+1]='.'; i+=2;
          for (k=n-lenFL/2+2; k<n; ++k, ++i) hstr[i]=hstr[k];
          hstr[i]=0;
      }
   }
   else { sprintf(hstr,"(null):%d", line); }

   for (i=strlen(hstr); i<lenFL; i++) hstr[i]=' ';

   strncpy(hstr+i,time_stamp+11,8);
   i+=8; hstr[i++]=' '; hstr[i]=0;

   if (tag && tag[0])
        sprintf(hstr+i, " %s ",tag);
   else { hstr[i]=' '; hstr[i+1]=0; }" ");
};


int check_update_header(
   char *hstr, const char* fmt,
   const char *istr=NULL, int lenFL=21
){
   int n,k=0,i=0;

   for (; fmt[i]; i++) {
      if (!isalpha(fmt[i]) && fmt[i]!='_' && fmt[i]!='.') {
         if (fmt[i]=='/') k=i+1; else break;
      }
   }

   if (fmt[i]!=':') return 0;
   n=i++;

   for (; fmt[i]; i++) if (!isdigit(fmt[i])) break;
   if (i==n+1 || isalpha(fmt[i])) return 0;

   if (istr && istr[0]) {
      strcpy(hstr,istr);
      n=strlen(hstr); hstr+=n; lenFL-=n;
   }

   if (i-k<lenFL) {
      memcpy(hstr,fmt+k,i-k);
      for (k=i-k; k<lenFL; k++) hstr[k]=' ';
   }
   else {
      memcpy(hstr,fmt+k,lenFL/2); k=lenFL/2;
      hstr[k]=hstr[k+1]='.';
      memcpy(hstr+k+2, fmt+i-k+2,k-2); k=2*k;
   }

   if (k<=lenFL) {
      for (; hstr[k]; k++) if (hstr[k]!=' ') break;
      if (k>lenFL) {
         for (n=k-lenFL; hstr[k]; k++) hstr[k-n]=hstr[k];
         hstr[k-n]=0;
      }
   }

   for (; fmt[i]; i++) {
      if (!isspace(fmt[i]) && fmt[i]!=':' && fmt[i]!='(' && fmt[i]!=')')
      break;
   }

   return i;
};


inline void wbSetLogLevel(unsigned l) {
    wblog(0, 0, "log_level", l);
}


int wblog(const char* file, int line, const char *fmt, ...) {

    static char log_header[128];
    static int log_level=1;
    static int then=0;

    char *s;

    time_t curr_time;
    struct tm *tblock;

    char fstr[128];      
    char time_stamp[32], tag[8];

    unsigned i,j,k,l,m,
         i2=0, isfmt=0, hflag=1, bflag=0, iflag=0, eflag=0;
    char c, *cp, cb;

    va_list ap;

    if (file==NULL && line<=0) {
       if (!strcmp(fmt,"log_level")) {
          va_start(ap,fmt);
          log_level=(unsigned)va_arg(ap,int);
          va_end(ap);
          return 0;
       }
    }

    if (log_level<0) return 0;

    if (fmt==NULL) wblog(FL,"ERR wblog requires input string / fmt.");
    if (fmt[0]=='\\') { hflag=0; fmt++; bflag=1; }

    strcpy(tag,"ERR"); eflag=wblog_checktag(fmt,tag); if (!eflag) {
    strcpy(tag,"WRN"); iflag=wblog_checktag(fmt,tag); if (!iflag) {
    strcpy(tag,"TST"); iflag=wblog_checktag(fmt,tag); if (!iflag) {
    strcpy(tag,"<i>"); iflag=wblog_checktag(fmt,tag); if (!iflag) {
    tag[0]=0; }}}
   #ifdef DBSTOP
      else { eflag=iflag; iflag=0; }
   #endif
    }

    if (log_level==0)
    if (!eflag && !iflag) return 0;

    if (eflag || iflag) fmt+=(eflag+iflag);
    if (eflag) printf("\n");

    curr_time = time(NULL);
    tblock = localtime(&curr_time);
    strcpy (time_stamp, asctime(tblock));

    if (then && then!=tblock->tm_mday) {
       char dstr[32];
       strftime(dstr,31,"%d-%b-%Y %T", tblock);
       printf("\n>> TODAY %s\n\n",dstr);



    }
    then=tblock->tm_mday;

    init_header(log_header,file,line,time_stamp,tag);


    for (;;) {
        if (fmt[0]!='%') break;
        if (fmt[1]!='N') break;

        printf("\n");
        fmt+=2;
    }

    if (hflag) { printf("%s",log_header); hflag=0; fflush(0); }

    i=isfmt=0;
    va_start(ap,fmt);

    while ((c=(*fmt++)))  { fstr[i++]=c; i2++;

        if (c=='%')  {
            isfmt=1;

            if (i>1) {
                if (fstr[i-2]=='%' || fstr[i-2]=='\\') {
                    isfmt=0;
                }
                else continue;
            }
        }


        if (fstr[i-1]==10)  {
            fstr[i-1] = 0;
            printf("%s\n",fstr); i2=0;
            fmt+=check_update_header(log_header,fmt);
            printf("%s",log_header); fflush(0);
            i=isfmt=0;
        }

        if (isfmt && c!='%')  {
            switch(c) {

                case '-':
                case '+':
                case '.':
                case ' ':
                case '0':
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                case 'L':           
                case 'l': break;    

                case 's':           
                     fstr[i]=0;
                     s=va_arg(ap, char*);

                     for (k=0; s[k];) {
                        if (s[k]=='\n') {
                           if (i>1) {
                              fstr[i-2]=0; printf("%s",fstr);
                              strcpy(fstr,"%s"); i=2;
                           }
                           printf("\n%s",log_header); fflush(0);
                           s+=(k+1); k=0;
                        }
                        else if (s[k]=='\\' && s[k+1]=='n') {
                           if (i>1) {
                              fstr[i-2]=0; printf("%s",fstr);
                              strcpy(fstr,"%s"); i=2;
                           }
                           printf("\n%s",log_header); fflush(0);
                           s+=(k+2); k=0;
                        }
                        else if (isspace(s[k])) { ++k; }
                        else break;
                     }

                     j=strlen(s);
                     if (j && (s[j-1]=='\\')) { bflag=2;
                        char s2[j+1]; strcpy(s2,s); s2[j-1]=0;
                        printf(fstr,s2);
                     }
                     else {
                        if (i2<=3) {"%s" as first format specs
                        k=check_update_header(log_header,s,"> ");
                        if (k) { s+=k;
                           printf("\r%s",log_header); fflush(0);
                           for (; fmt[0]; fmt++) {
                              if (!isspace(fmt[0]) && fmt[0]!=':'
                              && fmt[0]!='(' && fmt[0]!=')') break;
                           }
                        }}
                        printf(fstr,s);
                     }

                     i=isfmt=0;
                     break;

                case 'd':           
                case 'x':           
                case 'X':           
                     fstr[i]=0;
                     m=0; if (i>1) if (fstr[i-2]=='l') m++;
                     if (m) printf(fstr, va_arg(ap, long));
                     else   printf(fstr, va_arg(ap, int ));

                     i=isfmt=0;
                     break;

                case 'c':           
                     fstr[i]=0;
                     printf(fstr, va_arg(ap, int ));
                     i=isfmt=0;
                     break;

                case 'e':           
                case 'E':
                case 'g':
                case 'f':
                     fstr[i]=0;
                     m=0; if (i>1) if (fstr[i-2]=='L') m++;
                     if (m) printf(fstr, va_arg(ap, long double));
                     else   printf(fstr, va_arg(ap, double));
                     i=isfmt=0;
                     break;

                case 'N':           
                     fstr[i-2]=0;
                     printf("%s\n",fstr);
                     i=isfmt=0;
                     break;

                case 'D':           

                     fstr[i-2]=0; { char s[32]; unsigned l=0;
                     strncpy(s+l,time_stamp+0,10); l+=10;
                     strncpy(s+l,", ",         2); l+=2;
                     strncpy(s+l,time_stamp+20,4); l+=4; s[l]=0;
                     
                     printf("%s%s",fstr,s); }

                     i=isfmt=0;
                     break;

                case 'R':           
                     k = i;         
                     fstr[--i]=0;   

                     for (; int(i)>=0; i--)
                     if (fstr[i]=='%') break;

                     if (int(i)<0)
                     printf("%s:%d ERR Invalid format >%s<\n", FLINE, fmt-k);
                     else {
                        fstr[i] = 0;
                        printf("%s",fstr);

                        sscanf(fstr+i+1, "%d", &k);

                        sprintf(fstr, "%s", va_arg(ap, char *));

                        for (; k>0; k--)
                        printf("%s", fstr);
                     }

                     i=isfmt=0;
                     break;

                case 'B': 

                     k = i;
                     fstr[--i]=0;   

                     for (; int(i)>=0; i--)
                     if (fstr[i]=='%') break;

                     if (int(i)<0)
                     printf("%s:%d ERR Invalid format >%s<", FLINE, fmt-k);
                     else {
                        fstr[i] = 0;
                        printf("%s",fstr);    

                        i = sscanf(fstr+i+1, "%d", &k);
                        if (i<=0) k=8;    

                        m = va_arg(ap, int);    
                        cp = (char*)&m;         

                        if (k>0)
                        for (i=l=0; i<((k-1)/8+1); i++) {
                            for (cb=1, j=0; j<CHAR_BIT; ++j, cb<<=1) {
                                if (cp[i] & cb) printf ("1");
                                else            printf ("0");
                                if ((++l)>=k) break;
                            }
                        }
                     }

                     i=isfmt=0;
                     break;

                default:
                     fstr[i]=0;
                     strcat(fstr, " <ERR-F> ");
                     printf("%s", fstr);
                     fflush(stdout);
                     i=isfmt=0;
            }
        }
    }

    va_end(ap);

    fstr[i]=0;


    if (i && (fstr[i-1]=='\\')) { fstr[i-1]=0; bflag=2; }
    if (bflag<2) strcat(fstr,"\n");

    printf("%s",fstr);

    if (bflag || eflag || iflag) { fflush(0); doflush(); }
    if (eflag) { 
       ExitMsg();
    }

    return 1;
};


inline unsigned wblog_checktag(const char *istr, const char *tag) {

    const char *i,*j;
    unsigned n=strlen(tag);

    if (!n) return 0;
    i=strstr(istr,tag); if (i==0) return 0;

    for (j=i; j>istr; j--) if (isalnum(*j)) return 0;
    if (isalnum(i[n])) return 0;

    return unsigned((i+n)-j + (i[n] ? 1:0));
}


char wblog_findtoken(const char *istr, const char *tok, int maxoffset) {

    char f=1; const char *c=strstr(istr,tok);

    if (c==0) return 0;
    if (maxoffset>=0) if (c>istr+maxoffset) return 0;

    if (c>istr)
    if (isalnum(*(c-1))) f=0;
    if (isalnum(*(c+3))) f=0;

    return f;
}


#endif

