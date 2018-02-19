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

char USAGE[]="\n\
   msgsend <pid> <msg> [number];                        \n\
                                                        \n\
      <pid>  is message key of message queue to use     \n\
      <msg>  is text message to be sent                 \n\
      number optional value to be included with message \n\
                                                        \n\
   Wb,Jan20,08                                          \n\
";

#include "wblib.h"


void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    int i,pid=0; double val=0;
    wbstring msg;

    if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }
    if (nargin<2 || nargin>3 || nargout>2) { usage(); return; }

    mxGetNumber(argin[0],pid);
    msg.init(argin[1]);

    if (nargin>2)
    mxGetNumber(argin[2],val);

    i=WbMSGSnd(FL,pid,msg.data,val);

    if (nargout) argout[0]=numtoMx(i);
}


