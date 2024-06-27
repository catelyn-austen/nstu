#define _WINSOCK_DEPRECATED_NO_WARNINGS
#pragma comment(lib, "Ws2_32.lib")
#include <iostream>
#include <winsock2.h>
#include <thread>
#include <windows.h>
using namespace std;

void GetCMessage(SOCKET clientSock)
{
    char buff[512]{ };
    while (true)
    {
        int err = recv(clientSock, buff, 512, 0);
        if (err > 0)
        {
            for (int i = 0; buff[i] != '\0'; i++)
                cout << buff[i];
            cout << "> " << buff + 20 << endl;
        }
    }
}

int main()
{
    setlocale(LC_ALL, "ru");
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    WSADATA wsaData;
    WORD ver = MAKEWORD(2, 2);
    SOCKET clientSock;
    sockaddr_in servInfo;
    int err = 0;
    char buff[490]{ };


    err = WSAStartup(ver, &wsaData);
    if (err == 1)
    {
        cout << "WSAStartup failed." << endl;
        return 1;
    }

    clientSock = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP);
    if (clientSock == SOCKET_ERROR)
    {
        cout << "Unable to create socket" << endl;
        WSACleanup();
        return 1;
    }

    servInfo.sin_family = PF_INET;
    servInfo.sin_addr.s_addr = inet_addr("127.0.0.1");
    servInfo.sin_port = htons(2009);

    err = connect(clientSock, (sockaddr*)&servInfo, sizeof(servInfo));
    if (err == SOCKET_ERROR)
    {
        cout << "Unable to connect." << endl;
        WSACleanup();
        closesocket(clientSock);
        return 1;
    }
    cout << "Connected." << endl;

    while (true)
    {
        thread th(GetCMessage, clientSock);
        th.detach();
        cin.getline(buff, 490, '\n');
        err = send(clientSock, buff, 490, 0);
        if (err == SOCKET_ERROR)
            cout << "Unable to send message." << endl;

    }

    closesocket(clientSock);
    WSACleanup();
}