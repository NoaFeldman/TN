/* ------------------------------------------------------------------ //
   msgcheck.c -- receiving the above message
   (see msg_send.c routine for further information!)

   The essential points to note here are:

   The Message queue is opened with msgget (message flag 0666)
   and the same key as message_send.c.
   A message of the same type 1 is received from the queue
   with the message ``Did you get this?'' stored in rbuf.text.

   The full code listing for message_send.c's companion process,
   message_rec.c is as follows:

   Wb,Nov29,03

char USAGE[]="    \n\
   i=msgcheck();  \n\
   Wb,Jan18,08    \n\
";

#include "wblib.h"

wbMsg M(FL,"MEX");


void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    mxArray *S;

    wbMsgBuf rbuf;

    if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }

    M.Rec(FL,rbuf           );

    S=mxCreateStructMatrix(1,1,0,NULL);

    if (rbuf.time) {
       sprintf(str,"%s", ctime(&rbuf.time));
       str[strlen(str)-1] = 0;
    }
    else str[0]=0;

    mxAddField2Scalar(FL,S,"time", wbstring(str).toMx());
    mxAddField2Scalar(FL,S,"msg",  wbstring(rbuf.text).toMx());
    mxAddField2Scalar(FL,S,"val",  numtoMx(rbuf.val));
    mxAddField2Scalar(FL,S,"type", numtoMx(rbuf.type));

    mxAddField2Scalar(FL,S,"flag", numtoMx(
       (rbuf.type!=0 || rbuf.text[0]) ? 1 : 0
    ));


    argout[0]=S;
}

