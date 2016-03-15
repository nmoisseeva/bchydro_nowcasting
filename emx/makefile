CC = gcc
CFLAGS = 
IFLAGS  = -I/users/model/src/database/EMWX_lib \
		  -I/usr/lib/gcc/x86_64-redhat-linux/4.1.1/include/ \
          -I/usr/include/mysql
LFLAGS  = -L/users/model/src/database/EMWX_lib \
          -L/usr/lib64/mysql
LIBS   =  -lemwx_lib -lmysqlclient 

getstations_da.exe: getstations_da.o makefile
	$(CC) $(CFLAGS) getstations_da.o $(LFLAGS) $(LIBS) -o $@

getstations_da.o: getstations_da.c
	$(CC) $(CFLAGS) $(IFLAGS) -c getstations_da.c
