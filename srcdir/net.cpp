/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2014 Robert Baron, Charles P. Kollar,
  Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
  Daniel E. Weeks, and University of Pittsburgh

  This file is part of the Mega2 program, which is free software; you
  can redistribute it and/or modify it under the terms of the GNU
  General Public License as published by the Free Software Foundation;
  either version 3 of the License, or (at your option) any later
  version.

  Mega2 is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

  For further information contact:
      Daniel E. Weeks
      e-mail: weeks@pitt.edu

===========================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#if defined(_WIN) || defined(MINGW)
#include <winsock2.h>

#ifdef _WIN
#pragma comment (lib, "Ws2_32.lib")
#pragma comment (lib, "Mswsock.lib")
#pragma comment (lib, "AdvApi32.lib")
#endif /*  _WIN */

#include "cw_routines_ext.h"
#define SOCK_ERRh(sfd) (sfd == INVALID_SOCKET)
#define SOCK_ERR(sfd)  (sfd == SOCKET_ERROR)
#define ERR_STR() { errno = WSAGetLastError(); (errstr = "windows ..."); }
#define STR_ERR() /*LOCALFREE(errstr)*/
#define SOCK SOCKET
#define fdopen _fdopen
int errno;
const char *errstr;

#else /*  defined(_WIN) || defined(MINGW) */

#include <unistd.h>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <fcntl.h>
#define SOCK int
#define INVALID_SOCKET -1
#define SOCK_ERRh(sfd) (sfd < 0)
#define SOCK_ERR(sfd)  (sfd < 0)

// The system error message strings can be accessed directly using the external array 'sys_errlist'.
// The external value 'sys_nerr' contains a count of the messages in 'sys_errlist'.
// THE USE OF THESE VARIABLES IS DEPRECATED; strerror() or strerror_r() should be used instead.
// This fails for gcc version 4.5.2; target i386-pc-solaris2.11
#ifdef __sun
#define ERR_STR() errstr = strerror(errno)
#else /* __sun */
#define ERR_STR() errstr = errno < sys_nerr ? sys_errlist[errno] : "oops"
#endif /* __sun */

#define STR_ERR() 
const char *errstr;

#endif /*  defined(_WIN) || defined(MINGW) */

#define NOBOOLEAN
#include "common.h"

#include "error_messages_ext.h"
/*
     error_messages_ext.h:  errorf mssgf warnf
*/


void           test_socket_fd(void);
FILE          *http_response_body(const char *request, const char *host, unsigned short port);
void           socket_close_fd(FILE *FD);

void           test_socket_s(void);
typedef        void(*eat_fn_t)(int cnt, char *line, int len, int type, void *ap);
void           http_response(const char *request, const char *host, unsigned short port, eat_fn_t fn, void *ap);
void           show_response_line(int cnt, char *line, int len, int type);

void           socket_close_s(SOCK fd);
static int     refill(SOCK fd, char *buffer, size_t blen);

SOCK           http_request(const char *request, const char *host, unsigned short port);
SOCK           socket_fd(const char *host, unsigned short port);

static         char socket_buffer[4096];
//xx static         char socket_buffer[96];

struct shar {int line;} data;

void show_response(int cnt, char *line, int len, int type, void *ap)
{
    struct shar *sp = (struct shar *)ap;
    sp->line++;
    if (type > 0) msgvf("#%d, %s\n", sp->line, line);
}

void test_socket_s(void)
{
    data.line = 0;
    http_response("GET /pub/mega2/.version", "watson.hgen.pitt.edu", 80, show_response, &data);
}

static char *refill_overflow;
static char *refill_prev_end;

void http_response(const char *request, const char *host, unsigned short port, eat_fn_t fn, void *ap)
{
    SOCK fd;
    int status = 0;
    int len = 0;
    int type = 0;
    int cnt = 1;
    size_t line;
    char *linec, *linep;

    fd = http_request(request, host, port);
    if (SOCK_ERRh(fd)) return;

    refill_prev_end = NULL;
    linec = socket_buffer;
    linec[0] = 0;
    while (len >= 0) {
        if (len == 0) {
            len = refill(fd, socket_buffer, sizeof(socket_buffer));
            if (len < 0) break;
            linec = socket_buffer;
        }
        linep = strsep(&linec, "\n");
        len--;
        line = linec - linep - 1; /* -1 is for \n */
        len -= line;
        if (linec[-2] == '\r') {
            linec[-2] = 0;
            line--;
        }
        if (status == 0) {
            if (sscanf(socket_buffer, "HTTP/1.1 %3d", &status) != 1 ||
                status != 200) {
                if (status != 404 || debug)
                    warnvf("HTTP server returned %d %s\n", status, socket_buffer+12);
                socket_close_s(fd);
                return;
            }
        }
        if (debug) {
            msgvf("#%d: %s (L%d) (R%d) %d\n", cnt, linep, line, len, type);
        }

        fn(cnt, linep, (int)line, type, ap);

        if (*linep == 0) type++;
        cnt++;
    }
    return;
}

static int refill(SOCK fd, char *buffer, size_t blen)
{
    size_t len, overflow_size, total_length;

/*  blen=100; */        /* for testing */

    if (refill_prev_end == NULL) {
        refill_overflow = &buffer[blen];
        refill_prev_end = &buffer[blen];
    }
    if (refill_overflow < refill_prev_end) {
        overflow_size = refill_prev_end - refill_overflow;
        strncpy(buffer, refill_overflow, overflow_size);
    } else {
        overflow_size = 0;
    }

    len = recv(fd, buffer+overflow_size, blen-overflow_size, 0);
    if (len > 0) {
    } else if (len == 0) {
        buffer[0] = 0;
        socket_close_s(fd);
        return -1;
    } else {
        ERR_STR();
        switch(errno) {
        default:
            warnvf("test: recv() %d failed with errno %d (\"%s\")\n",
		   len, errno, errstr);
            STR_ERR();
            socket_close_s(fd);
            return -1;
        }
    }

    total_length = overflow_size + len;
    refill_prev_end = buffer + total_length;
    if (total_length == blen) {
        refill_overflow = buffer + total_length;
        while (--refill_overflow > buffer && *refill_overflow != '\n') ;
        refill_overflow++;
    } else {
        refill_overflow = &buffer[total_length];
        buffer[total_length] = 0;
    }
    overflow_size = refill_prev_end - refill_overflow;
    return (int)total_length - overflow_size;
}

void test_socket_fd(void)
{
    FILE *FD = http_response_body("GET /pub/mega2/.version", "watson.hgen.pitt.edu", 80);
    size_t len;
    char *ret;

    if (FD == NULL) return;
    while (!feof(FD)) {
        ret = fgets(socket_buffer, sizeof(socket_buffer), FD);
        if (ret == NULL) {
            if (feof(FD)) {
                socket_close_fd(FD);
                return;
            }
            ERR_STR();
            switch(errno) {
            default:
                warnvf("test: fgets() %s failed with errno %d (\"%s\")\n",
		       ret, errno, errstr);
                STR_ERR();
                socket_close_fd(FD);
                return;
            }
        } else {
            if (debug) {
                len = strlen(socket_buffer);
                socket_buffer[len-1] = 0;
                mssgf(socket_buffer);
            }
        }
    }
}

FILE *http_response_body(const char *request, const char *host, unsigned short port)
{
    FILE *FD;
    int status = 0;
    size_t len;
    char *ret;

    int fd = http_request(request, host, port);
    if (fd < 0) return NULL;
    FD = fdopen(fd, "r");
    if (FD == NULL) {
        return FD;
    }

    while (!feof(FD)) {
        ret = fgets(socket_buffer, sizeof(socket_buffer), FD);
        if (ret == NULL) {
            if (feof(FD)) {
                socket_close_fd(FD);
                return NULL;
            }
            ERR_STR();
            switch(errno) {
            default:
                warnvf("test: fgets() %s failed with errno %d (\"%s\")\n",
		       ret, errno, errstr);
                STR_ERR();
                socket_close_fd(FD);
                return NULL;
            }
        } else {
            if (status == 0) {
                if (sscanf(socket_buffer, "HTTP/1.1 %3d", &status) != 1 ||
                    status != 200) {
                    if (status != 404 || debug)
                        warnvf("HTTP server returne d %d %s\n", status, socket_buffer+12);
                    socket_close_fd(FD);
                    return NULL;
                }
            }
            len = strlen(socket_buffer);
            if (debug) {
                socket_buffer[len-1] = 0;
                mssgf(socket_buffer);
            }
            if (len == 2) return FD;
        }
    }
    return NULL;
}

SOCK http_request(const char *request, const char *host, unsigned short port)
{
    const char *req = "%s HTTP/1.1\r\nHost: %s\r\nConnection: close\r\n\r\n";
    long len;

    SOCK fd = socket_fd(host, port);
    if (SOCK_ERRh(fd)) return INVALID_SOCKET;

#if defined(_WIN) || defined(MINGW)
    int soptbuff = 5 * 1000;
#else
    struct timeval soptbuff = { 5, 0};
#endif
    int sret;
    sret =           setsockopt(fd, SOL_SOCKET, SO_SNDTIMEO, (const char *)&soptbuff, sizeof soptbuff);
    if (sret != 0 || setsockopt(fd, SOL_SOCKET, SO_RCVTIMEO, (const char *)&soptbuff, sizeof soptbuff) != 0 ) {
        ERR_STR();
        warnvf("setsocketopt: SND/RCV timeout failed: errno %d (%s)\n", errno, errstr);
    }

    sprintf(socket_buffer, req, request, host);
    len  = strlen(socket_buffer);

    len = send(fd, socket_buffer, len, 0);
    if (SOCK_ERR(len)) {
        ERR_STR();
        switch(errno) {
        default:
            warnvf("http_request: write() failed with errno %d (\"%s\")\n",
		   errno, errstr);
            STR_ERR();
            socket_close_s(fd);
            return INVALID_SOCKET;
        }
    }
    return fd;
}

SOCK socket_fd(const char *host, unsigned short port)
{
    int err;
    SOCK  sfd;
    struct sockaddr_in saddr;
    struct hostent *mega2_server;

#if defined(_WIN) || defined(MINGW)
    WSADATA wsaData;
    WORD    version = MAKEWORD( 2, 2 ); /* current 2.2 */

    err = WSAStartup( version, &wsaData );
    if (err != 0) {
        warnvf("socket_fd: WSAStartup() failed with return value %d\n", WSAGetLastError());
        return INVALID_SOCKET;
    }
#endif

    sfd = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
    if (SOCK_ERRh(sfd)) {
        ERR_STR();
        switch(errno) {
        default:
            warnvf("socket_fd: socket() failed with errno %d (\"%s\")\n", errno, errstr);
            STR_ERR();
#if defined(_WIN) || defined(MINGW)
            WSACleanup();
#endif
            return INVALID_SOCKET;
        }
    }

    mega2_server = gethostbyname(host);
    if (mega2_server == NULL) {
        warnvf("socket_fd: gethostbyname(%s) failed\n", host);
        socket_close_s(sfd);
        return INVALID_SOCKET;
    }

    memset(&saddr, 0, sizeof (struct sockaddr_in));
    saddr.sin_family      = AF_INET;
    saddr.sin_port        = htons(port);

    memcpy(&saddr.sin_addr, mega2_server->h_addr_list[0], sizeof (struct in_addr));

#if defined(_WIN) || defined(MINGW)
    int soptbuff = 1;
    if (ioctlsocket(sfd, FIONBIO, (u_long *)&soptbuff) != NO_ERROR) {
        ERR_STR();
        warnvf("ioctlsocket: set non blocking failed\n");
    }
#else
    int fl = fcntl(sfd, F_GETFL);
    fcntl(sfd, F_SETFL, O_NONBLOCK);
#endif

    err = connect(sfd, (struct sockaddr *)&saddr, sizeof (struct sockaddr_in));

    fd_set fds;
    FD_ZERO(&fds);
    FD_SET(sfd, &fds);

    struct timeval sopttm = { 5, 0 };
    err = select(sfd + 1, NULL, &fds, NULL, &sopttm);

    if (SOCK_ERR(err) || err == 0) {
#if defined(_WIN) || defined(MINGW)
        if (err == 0) {
            errno = 0;
            errstr = "Operation timed out";
        } else
#else
        if (err == 0) {
            errno = 60;
            err   = -1;
        } // fall thru
#endif
            ERR_STR();
        switch(errno) {
        default:
            fflush(stdout); fflush(stderr);
            warnvf("socket_fd: connect(%s:%d) failed with errno %d (\"%s\")\n",
		   host, port, errno, errstr);
            STR_ERR();
            socket_close_s(sfd);
            return INVALID_SOCKET;
        }
    }

#if defined(_WIN) || defined(MINGW)
    soptbuff = 0;
    if (ioctlsocket(sfd, FIONBIO, (u_long *)&soptbuff) != NO_ERROR) {
        ERR_STR();
        warnvf("ioctlsocket: set non blocking failed\n");
    }
#else
    fcntl(sfd, F_SETFL, fl);
#endif
    return sfd;
}

void socket_close_s(SOCK fd)
{
#if defined(_WIN) || defined(MINGW)
    closesocket(fd);
    WSACleanup();
#else
    close(fd);
#endif
}

void socket_close_fd(FILE *FD)
{
#if defined(_WIN) || defined(MINGW)
    int fd = FD->_file;
    closesocket(fd);
    WSACleanup();
#else
    fclose(FD);
#endif
}

/* FILE *http_request(const char *request, const char *host, int port) */
/* { */
/*     char *req = "%s HTTP/1.1\r\nHost: %s\r\nConnection: close\r\n\r\n"; */
/*     int len; */

/*     int fd = socket_fd(host, port); */
/*     if (fd < 0) { */
/* #ifdef _WIN */
/*         WSACleanup(); */
/* #endif */
/*         return NULL; */
/*     } */

/*     sprintf(socket_buffer, req, request, host); */
/*     len  = strlen(socket_buffer); */

/*     len = write(fd, socket_buffer, len); */
/*     if (len < 0) { */
/*         switch(errno) { */
/*         default: */
/*             errorvf("http_request: write() failed with errno %d (\"%s\")\n", */
/*                     errno, errno < sys_nerr ? sys_errlist[errno] : "oops"); */

/* #ifdef _WIN */
/*             closesocket(fd); */
/*             WSACleanup(); */
/* #else */
/*             close(fd); */
/* #endif */
/*             return NULL; */
/*         } */
/*     } */

/*     return fdopen(fd, "r"); */
/* } */

/* int socket_fd(const char *host, int port) */
/* { */
/*     int err; */
/*     int sfd; */
/*     struct sockaddr_in saddr; */
/*     struct hostent *mega2_server; */

/* #ifdef _WIN */
/*     WSADATA wsaData; */
/*     WORD version = MAKEWORD( 2, 0 ); */

/*     err = WSAStartup( version, &wsaData ); */
/* #endif */

/*     sfd = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP); */
/*     if (sfd < 0) { */
/*         switch(errno) { */
/*         default: */
/*             errorvf("socket_fd: socket() failed with errno %d (\"%s\")\n", */
/*                     errno, errno < sys_nerr ? sys_errlist[errno] : "oops"); */
/*             return -1; */
/*         } */
/*     } */

/*     mega2_server = gethostbyname(host); */
/*     if (mega2_server == NULL) { */
/*         errorvf("socket_fd: gethostbyname(%s) failed\n", host); */

/*         close(sfd); */
/*         return -1; */
/*     } */

/*     memset(&saddr, 0, sizeof (struct sockaddr_in)); */
/*     saddr.sin_family      = AF_INET; */
/*     saddr.sin_port        = htons(port); */

/*     memcpy(&saddr.sin_addr, mega2_server->h_addr_list[0], sizeof (struct in_addr)); */

/*     err = connect(sfd, (struct sockaddr *)&saddr, sizeof (struct sockaddr_in)); */
/*     if (err < 0) { */
/*         switch(errno) { */
/*         default: */
/*             errorvf("socket_fd: connect(%s:%d) failed with errno %d (\"%s\")\n", */
/*                     host, port, errno, errno < sys_nerr ? sys_errlist[errno] : "oops"); */

/*             close(sfd); */
/*             return -1; */
/*         } */
/*     } */

/*     return sfd; */
/* } */
