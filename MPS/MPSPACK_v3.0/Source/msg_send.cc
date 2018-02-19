/* ------------------------------------------------------------------ //
   Example: Sending messages between two processes
   http://www.cs.cf.ac.uk/Dave/C/node25.html#SECTION002540000000000000000

   The following two programs should be compiled and run at the same time
   to illustrate basic principle of message passing:

   message_send.c
   -- Creates a message queue and sends one message to the queue.

   message_rec.c
   -- Reads the message from the queue.

   Wb,Nov29,03






#define MAIN
#define WBLOG_NO_INIT_MSG

char USAGE[256];

#include <wblib.h>

int main (int argc, char **argv)  {

    sprintf (USAGE,
       "\nUSAGE: %s <pid> <variable|string> [number]\n\n",
        argv[0]
    );
    if (argc!=3 && argc!=4) { usage(); return -1; }

    float fl=0.;
    int i, pid=0;


    i=sscanf(argv[1],"%d", &pid);
    if (i!=1) { usage(); return -1; }

    if (argc>3) {
        i=sscanf(argv[3],"%f", &fl);
        if (i!=1) { usage(); return -1; }
    }


    WbMSGSnd (FL, pid, argv[2], fl);

    return 0;
}

