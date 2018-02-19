#ifndef __WB_LOG_HH__
#define __WB_LOG_HH__

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// class WbLog that enforces proper unsetting of WBLOG_ERR_CONT
// after local environment terminates // Wb,Feb24,11

class WbLog {

  public:

    WbLog(const char *F, int L, const char *s)
     : err_cont(0) {
       if (!s || !s[0]) return;

       if (!strcmp(s,"continue on ERR")) {
          WBLOG_ERR_CONT=1;
          if (WBLOG_ERR_COUNT) { wblog(F_L,"WRN %s() "
             "got wblog error count %d !??",FCT,WBLOG_ERR_COUNT); 
             WBLOG_ERR_COUNT=0;
          }
       }
       else wblog(F_L,"ERR %s() invalid setting '%s'",FCT,s);
    };

   ~WbLog() {
       check(); WBLOG_ERR_CONT=0; 
    }

    void check(const char *F=0, int L=0) {
       if (!WBLOG_ERR_COUNT) { return; }

       unsigned e=WBLOG_ERR_COUNT;

       if (err_cont) {
          if (e==1) wblog(F_L,"ERR"); else
          if (e) wblog(F_L,"ERR %d errors",e);
       }
       else {
          if (e==1) wblog(F_L,"WRN"); else
          if (e) wblog(F_L,"WRN %d errors",e);
       }
    };

    bool ok() {
       if (WBLOG_ERR_COUNT) return 0;
       return 1;
    }

  protected:
  private:

    char err_cont;
};


#endif

