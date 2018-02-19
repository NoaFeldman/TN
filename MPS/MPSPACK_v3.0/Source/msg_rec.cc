/* ------------------------------------------------------------------ //
   message_rec.c -- receiving the above message
   (see msg_send.c routine for further information!)

   The essential points to note here are:

   The Message queue is opened with msgget (message flag 0666)
   and the same key as message_send.c.
   A message of the same type 1 is received from the queue
   with the message ``Did you get this?'' stored in rbuf.mtext.

   The full code listing for message_send.c's companion process,
   message_rec.c is as follows:

   Wb,Nov29,03

char USAGE[]="\n\
  Waiting for messages to this process using its PID        \n\
  as message queue ID. Send Message 'STOP' to terminate ... \n\
  \n\
";

#define MAIN
#define WBLOG_NO_INIT_MSG

#include <wblib.h>

int main() {

    wbMsg M(FL,"TST_Receive");
    wbMsgBuf rbuf;

    printf(USAGE);

    for (int i=0; ; i++) {
        M.Rec(FL,rbuf           );
        if (!strcmp(rbuf.text, "STOP")) break;

        sleep(1);
        printf("%3d: >%s< %d %g \r",i, rbuf.text, (int)rbuf.type, rbuf.val);
        fflush(0);
    }

    return 0;
}

