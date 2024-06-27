using System.IO;

public class Program
{
    public static string Shuffle()
    {
        char[] alphabet = "абвгдеёжзийклмнопрстуфхцчшщъыьэюя".ToCharArray();
        Random r = new Random();
        int j;
        for (int i = alphabet.Length - 1; i > 0; i--)
        {
            j = r.Next(i + 1);
            (alphabet[i], alphabet[j]) = (alphabet[j], alphabet[i]);
        }
        return new string(alphabet);
    }
    public static void Main()
    {
        string alphabet = "ктёдхяншэгибвйоюжчаурсфъмлпзьыцещ";
        string key = "абвгде_жзийклмнопрстуфх_чшщ_ыьэюя";
        string essay;
        StreamWriter sw = new StreamWriter("key.txt");
        sw.WriteLine(key);
        sw.Close();

        StreamReader sr = new StreamReader("in.txt"); // в файле эссе
        sw = new StreamWriter("out.txt"); // файл на выход

        essay = sr.ReadToEnd();
        sw.Write("{\n" + "\"essay\": " + "\"" + essay +"\",\n" +
                         "\"key\": " + "\"" + key + "\",\n" +
                         "\"encryptedEssay\": " + "\"");
        essay = essay.ToLower(); // меняем регистр
        foreach (char c in essay) // удаляем лишние символы
        {
            if (!alphabet.Contains(c))
                essay = essay.Replace(c.ToString(), string.Empty);
        }
        
        foreach (char c in essay) // шифруем
        {
            sw.Write(key[alphabet.IndexOf(c)]);
        }
        sw.Write("\"\n" + "}");
        sr.Close();
        sw.Close();
    }
}