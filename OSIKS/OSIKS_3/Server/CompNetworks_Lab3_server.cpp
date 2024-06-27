#define _WINSOCK_DEPRECATED_NO_WARNINGS
#pragma comment(lib, "Ws2_32.lib")
#include <iostream>
#include <winsock2.h>
#include <omp.h>
#include <vector>
#include <thread>
#include <windows.h>
using namespace std;



void ClientThread(SOCKET ClientSocket, sockaddr_in from, vector<SOCKET>& ClientSocketArr)
{
    char buff[512]{};
    string ip = inet_ntoa(from.sin_addr);
    string s;
    while (true)
    {
        for (int i = 0; i < ip.size(); i++)
            buff[i] = ip[i];
        int err = recv(ClientSocket, &buff[20], 490, 0);
        if (err > 0)
        {
            s = &buff[20];
            if (s == "end")
                break;
            cout << ip << " > " << &buff[20] << endl;
        }
        else
        {
            int i = 0;
            for (; i < ClientSocketArr.size(); i++)
                if (ClientSocketArr[i] == ClientSocket)
                    break;
            ClientSocketArr.erase(ClientSocketArr.begin() + i);
            break;
        }
        for (int i = 0; i < ClientSocketArr.size(); i++)
        {
            if (ClientSocketArr[i] != ClientSocket)
                send(ClientSocketArr[i], buff, 512, 0);
        }
    }
    cout << ip << " disconnected" << endl;
    closesocket(ClientSocket);
}



int main()
{
    setlocale(LC_ALL, "ru");
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    WSADATA wsaData;
    WORD ver = MAKEWORD(2, 2);
    SOCKET ListenSock, ClientSock;
    vector<SOCKET> ClientSockArr;
    sockaddr_in servAddr;
    int err = 0;
    char buff[512]{ };


    err = WSAStartup(ver, &wsaData);
    if (err == 1)
    {
        cout << "WSAStartup failed." << endl;
        return 1;
    }

    ListenSock = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP);
    if (ListenSock == INVALID_SOCKET)
    {
        cout << "Unable to create server socket" << endl;
        WSACleanup();
        return SOCKET_ERROR;
    }
    servAddr.sin_family = PF_INET;
    servAddr.sin_addr.s_addr = inet_addr("127.0.0.1");
    servAddr.sin_port = htons(2009);

    err = bind(ListenSock, (sockaddr*)&servAddr, sizeof(servAddr));
    if (err == SOCKET_ERROR)
    {
        cout << "Bind socket error." << endl;
        closesocket(ListenSock);
        WSACleanup();
        return 1;
    }

    err = listen(ListenSock, 50);
    if (err == SOCKET_ERROR)
    {
        cout << "Listen socket error." << endl;
        closesocket(ListenSock);
        WSACleanup();
        return 1;
    }

    cout << "Server started." << endl;
    while (true)
    {
        sockaddr_in from;
        int fromlen = sizeof(from);
        ClientSock = accept(ListenSock, (sockaddr*)&from, &fromlen);
        if (ClientSock != INVALID_SOCKET)
        {
            cout << "IP " << inet_ntoa(from.sin_addr) << " connected" << endl;
            ClientSockArr.push_back(ClientSock);
        }

        thread th(ClientThread, ClientSock, from, ref(ClientSockArr));
        th.detach();
    }
}