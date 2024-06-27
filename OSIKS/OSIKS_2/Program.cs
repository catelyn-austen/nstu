using System;
using System.IO;
using System.Net;
using System.Net.Sockets;
using System.Text;

class TCP_Server
{
    public static void Main()
    {
        TcpListener server = null;
        try
        {
            Int32 port = 2009; // порт сервера
            IPAddress localAddr = IPAddress.Parse("25.37.214.127"); // ip-адрес сервера (интерфейс)
            server = new TcpListener(localAddr, port); //TcpListener - класс TCP-сервера из .Net Framework Class Library
            server.Start(); // начинаем ожидание подсоединений клиентов на интерфейсе localAddr и порту port
            Byte[] bytes = new Byte[1000]; // буффер для приема сообщений
            String data; // строка для вывода на экран
            String answer_message; // ответ клиенту
            
            while (true) // цикл обработки подсоединений клиентов
            {
                Console.Write("Waiting for a connection... ");
                TcpClient client = server.AcceptTcpClient(); // ждем клиента и его присоединения
                Console.WriteLine("Connected!");
                NetworkStream stream = client.GetStream(); // вводим поток stream для чтения и записи через установленное соединение
                int i = stream.Read(bytes, 0, bytes.Length); // как и в программе клиента stream.Read возвращает значение 0,
                                                             // если запрошено 0 байтов или если в потоке больше нет данных
                if (i > 0) // соответственно, если данные есть
                {
                    data = System.Text.Encoding.ASCII.GetString(bytes, 0, i); // преобразуем принятые данные в строку ASCII string.
                    Console.WriteLine("Received: {0}", data);
                    int Count = data.Length; // считаем количество символов в строке
                    answer_message = data + " (" + Count.ToString() + ")";
                    Console.WriteLine("Sent: {0}", answer_message);
                    byte[] msg = System.Text.Encoding.ASCII.GetBytes(answer_message); // преобразуем строчку-ответ сервера в массив байт
                    stream.Write(msg, 0, msg.Length); // отправляем ответ
                }
                client.Close();
            }
        }
        catch (SocketException expt)
        {
            Console.WriteLine("SocketException: {0}", expt);
        }
        finally
        {
            server.Stop(); // Stop listening for new clients.
        }
        Console.WriteLine("\nHit enter to continue...");
        Console.Read();
    }
}