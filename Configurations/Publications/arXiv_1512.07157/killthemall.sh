#!/bin/bash

USERNAME=paulayan

#to kill all the jobs
qstat -u $USERNAME | grep "$USERNAME" | cut -d"." -f1 | xargs qdel

#to kill all the running jobs
#qstat -u $USERNAME | grep "R" | cut -d"." -f1 | xargs qdel

#to kill all the queued jobs
#qstat -u $USERNAME | grep "Q" | cut -d"." -f1 | xargs qdel
