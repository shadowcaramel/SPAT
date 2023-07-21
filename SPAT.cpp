#define _USE_MATH_DEFINES // для M_PI - числа пи
#include <iostream>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <iomanip> // для setw...
#include <string>  // подключаем строки
#include <fstream> // подключаем файлы
#include <stdlib.h>
//#include <math.h>
#include <cmath>
#include <complex> // для дифр поправки для коэф. затухания
#include <vector> // для функции автокорреляции (пока лишь)
#include <omp.h> // для быстрого подсчета автокорреляции

#include <stdio.h> 
#include <time.h> 


using namespace std;

struct dots // элемент данных
{
    double t; // время
    double A; // амплитуда
};

struct dwtresult // вейвлет-преобразование в точке W(a,b)
{
    double a = 0;          // масштаб
    double b = 0;          // временной сдвиг
    complex<double> W; // результат вейвлет-преобразования 
};

vector<filesystem::path> get_filenames(filesystem::path directory_path) // получаем имена файлов с данными
{
    vector<filesystem::path> filenames;
    for (const auto& entry : filesystem::directory_iterator(directory_path)) // directory_path директория с интересующими файлами
    {
        //cout << entry.path() << endl; // выводим пути найденных файлов
        filenames.push_back(entry.path()); // загоняем в вектор с именами файлов
    }
    return filenames;
} 

vector<dots> fill(filesystem::path filename, bool ignore5) // заполнение вектора данных собственно данными из файла
{
    vector<dots> v; // сюда складываем значения
    ifstream file1(filename); // файл, из которого читаем
    string s; // сюда будем класть считанные строки

    if (!file1)
        cout << "error reading file" << endl;    

    if (ignore5 == true) // игнорируем  первые 5 строк, которые есть в изначальном файле с осциллограммой
    {
        file1.ignore((numeric_limits<streamsize>::max)(), '\n');
        file1.ignore((numeric_limits<streamsize>::max)(), '\n');
        file1.ignore((numeric_limits<streamsize>::max)(), '\n');
        file1.ignore((numeric_limits<streamsize>::max)(), '\n');
        file1.ignore((numeric_limits<streamsize>::max)(), '\n');
    }
    
    //считываем cтроки с файла и записываем в вектор
    while (getline(file1, s)) {                
        auto pos = s.find(",");
        if (pos != string::npos)
        {
            v.push_back({ stod(s.substr(0, pos)), stod(s.substr(pos + 1)) }); //string to double // время, амплитуда добавляется как точка            
        }                    
    }
    file1.close(); // закрываем файл

    return v;
    //считали, заполнили
}

unsigned __int64 findmax(vector<dots> v) // ищем глобальный максимум амплитуды
{   
    double max = v[0].A; // максимум амплитуды
    unsigned __int64 maxindex = 0;
    unsigned __int64 i; //номер макс элемента
    for (i = 1; i < v.size(); i++) {
        if (v[i].A == 0) continue; // есть ли смысл?
        if (v[i].A > max) {
            max = v[i].A;
            maxindex = i;
        }
    }
    return maxindex;
}

void del(vector<dots> &a, unsigned __int64 ind, unsigned __int64 deltai) // другой вариант (лучше, вроде) зануления всех элементов вокруг выбранного элемента в радиусе deltai (по индексам)
{
    // deltai - "радиус" в индексах для зануления элементов
    double coef = 2.0;
    unsigned __int64 upper_bound = ind + coef * deltai;
    unsigned __int64 lower_bound = ind - coef * deltai; // при unsigned int возникает ошибка, если получается отрицательное (ind-deltai) при использовании эквивалентного условия

    if (ind < deltai) // <=> (ind - deltai < 0)
        lower_bound = 0;

    unsigned __int64 maxind = a.size(); // максмальный индекс
    if (upper_bound > maxind)
        upper_bound = maxind; //чтобы не выйти за пределы
    
    for (unsigned __int64 i = lower_bound; i < upper_bound; i++)
        a[i].A = 0; // зануляем <=> удаляем
}

dots average(vector<dots> a) // подсчет среднеарифметического значения вектора
{
    unsigned __int64 maxind = a.size();
    dots avgvalue = { 0,0 };
    for (unsigned __int64 i = 0; i < maxind; i++)
    {
        avgvalue.t += a[i].t; // вычисляем среднее время (незачем)
        avgvalue.A += a[i].A; // и амплитуду
    }
    avgvalue.t = avgvalue.t / maxind;
    avgvalue.A = avgvalue.A / maxind;

    return avgvalue;
}

void sorta(vector<dots> a) { // сортировка массива пузырьком по возрастанию
    unsigned __int64 maxind = a.size();
    double temp0, temp1; //временные переменные для обмена, соотв. a[...].t и a[...].A
    for (unsigned __int64 i = 0; i < maxind - 1; i++) {
        for (unsigned __int64 j = 0; j < maxind - i - 1; j++) {
            if (a[j].t == 0 || a[j].A == 0) continue;
            if (a[j].t > a[j + 1].t) {
                temp0 = a[j].t;
                temp1 = a[j].A;

                a[j].t = a[j + 1].t;
                a[j].A = a[j + 1].A;

                a[j + 1].t = temp0;
                a[j + 1].A = temp1;
            }
        }
    }
}

void sortd(vector<dots> a) { // сортировка массива пузырьком по убыванию
    unsigned __int64 maxind = a.size();
    double temp0, temp1; //временные переменные для обмена, соотв. a[0][...] и a[1][...]
    for (unsigned __int64 i = 0; i < maxind - 1; i++) {
        for (unsigned __int64 j = 0; j < maxind - i - 1; j++) {
            if (a[j].t == 0 || a[j].A == 0) continue;
            if (a[j].t < a[j + 1].t) {
                temp0 = a[j].t;
                temp1 = a[j].A;

                a[j].t = a[j + 1].t;
                a[j].A = a[j + 1].A;

                a[j + 1].t = temp0;
                a[j + 1].A = temp1;
            }
        }
    }
}

double k_La(double f, double a, double c_L) {  // параметр k_La для расчета дифр. поправки; с_L - скорость; a - диаметр(?) электрода; f - частота
    return 2 * M_PI * f * a / c_L;
}

double iterative_avg(vector<dots> a, unsigned int k, double eps, unsigned int n) { //функция для подсчета "среднего" уровня шума
    // итеративно вычисляем "среднее" занчение, отбрасывая большие (в k раз больше среднего значения) значения

    // a - массив данных для обработки    
    // k - коэф. для сравнения среднего и текущего значений; используется для отфильтровывания больших значений, т.е. не шума
    // eps - число для оценки изменений нового и старого средних; "точность"
    // n - максимальное число итераций

    unsigned __int64 length = a.size();
    double avg = 0;
    for (unsigned __int64 i = 0; i < length; i++)
        avg += abs(a[i].A);
    avg = avg / length; // находим "среднее" арифм по абс. значением

    cout << "iterative avg: avg=" << avg << endl;

    double newavg = 0;

    for (unsigned __int64 i = 0; i < length; i++) { // считаем newavg для использования в условии ниже
        if (abs(a[i].A) > avg * k)
            newavg += avg * k;
        else newavg += abs(a[i].A);
    }
    newavg = newavg / length;

    cout << "iterative avg: newavg=" << newavg << endl;

    double ratio = avg / newavg; //для проверки условия

    cout << "iterative avg: initial avg/newavg=" << avg / newavg << endl;

    for (unsigned int i = 0; ((i < n) && (ratio > 1 + eps)); i++) {

        cout << "iterative avg: iteration #" << i << endl;

        newavg = 0;

        for (unsigned __int64 i = 0; i < length; i++) { // на основе ранее вычисленного среднего вычисляем новое "модифицированное" среднее
            if (abs(a[i].A) > avg * k)
                newavg += avg * k;
            else newavg += abs(a[i].A);
        }
        newavg = newavg / length;

        cout << "iterative avg: newavg=" << newavg << endl;

        ratio = avg / newavg;

        cout << "iterative avg: avg/newavg=" << ratio << endl;

        avg = newavg;


    }
    return newavg;
}

double autocorrelation(vector<dots> a, unsigned __int64 j) { // функция автокорреляции сигнала
    // a - массив данных {t,A}
    // j - сдвиг по времени (индексу)
    // n - число элементов в дискретном сигнале
    unsigned __int64 n = a.size();
    double R = 0;

    // значение по методу трапеций, опять же, очень близко к значению по методу прямоугольников
    for (unsigned __int64 i = 0; i < n - j - 1; i++) { // метод трапеций 
        R += (a[i].A * a[i + j].A + a[i + 1].A * a[i + j + 1].A) / 2 * (a[i + 1].t - a[i].t);
        //Q++; // считаем кол-во интегральных слагаемых
    }
    R = R / (a[n - 1].t - a[0].t);   

    //cout << "R=" << R << endl;
    return R;
}

double norma(vector<dots> array) { // норма функции сигнала u(t) для нормирования ф-и автокорреляции
    // вычисляем интеграл N=  1/{t=t(maxindex)} * int_{t=0}^{t=t(maxindex)} u(t)^2 dt методом трапеций
    unsigned __int64 maxindex = array.size();

    double N = 0;
    for (unsigned __int64 i = 0; i < maxindex - 1; i++) // -2, т.к. метод трапеций
        N += (array[i].A * array[i].A + array[i + 1].A * array[i + 1].A) / 2 * (array[i + 1].t - array[i].t);
    

    return N / (array[maxindex - 1].t - array[0].t);
}

vector<dots> autocorrelation_calc(vector<dots> array)  // вычисляем функцию автокорреляции сигнала, записываем в файл
{
    unsigned __int64 maxindex = array.size();
    cout << "autocorrelation_calc: maxindex = " << maxindex << endl;
    //double** acor = new double* [2];
    //acor[0] = new double[maxindex];
    //acor[1] = new double[maxindex];
    vector<dots> acor(maxindex);

    double u_norma = norma(array);

    cout << "u_norma=" << u_norma << endl;

    //cout << "autocor=" << autocorrelation(u, 0, maxindex);

#pragma omp parallel for schedule (dynamic, 100) //num_threads(8)    //dynamic, 1000
    //for ( unsigned __int64 i = 0; i < counter - 1; i++) {
    // записываем в массив acor результаты ф-и автокорреляции по аналогии с массивом u: acor[0][...] - время,  acor[1][...] - значение ф-и автокор.
    for (__int64 i = 0; i < static_cast<__int64>(maxindex); i++) 
    {
        // вычисляем функцию автокорреляции R(tau)=1/T int_0^T  u(t)*u(t+tau) dt
        acor[i].t = array[i].t - array[0].t;
        acor[i].A = autocorrelation(array, i) / u_norma;
        //acor.push_back({ array[i].t - array[0].t, autocorrelation(array, i) / u_norma });
        //cout << "acor=" << acor[1][i] << endl;        
    }    
    return acor;
}

unsigned __int64 delta_in_indices(double delta, double duration, unsigned __int64 maxindex) { // преобразование "ширины" импульса в ширину "по индексам", т.е. время delta занимает столько-то индексов
    // duration - длительность осциллограммы в секундах
    // maxindex - длина массива

    return static_cast<unsigned __int64> ((delta * maxindex) / duration);
}

// убрано
/*unsigned __int64 deltai_modification(unsigned __int64 deltai, double frequency) { // зависимость deltai (ширины импульса) от неких параметров, пока что только от частоты
    return static_cast<unsigned __int64> (ceil(deltai * 15 * pow(frequency * 1e-6, -0.8))); // чем меньше частота, тем "шире" пик, и для компенсации этого увеличиваем k на низких частотах
}*/

double k_modification(double k, double f)  // зависимость k (что-то типа отношения сигнал/шум) от внешних параметров, пока что от частоты
{
    return k *= -0.02 * f * f * 1e-12 + 1.15 * f * 1e-6 - 1; // ~ чем больше частота, тем "выше" пик    
}

void filter_noise(vector<dots> &a, double k, double avg)  // "фильтрация" массива амплитуд от "шума" (небольших амплитуд) на основе коэф. k и "среднего" значения        
{
    for (unsigned __int64 i = 0; i < a.size(); i++) 
    {
        //u[1][i] = abs(u[1][i]); // можно применить из-за несимметричности (иногда) осциллограммы относительно оси абсцисс
        /*a[i].A -= avg * k;
        if (a[i].A < 0) // если за вычетом avg*k получилось отриц. значение, то зануляем его, и в дальнейшем интересовать нулевые значения не будут
            a[i].A = 0;*/
        if ((a[i].A > -1*k*avg) && (a[i].A < k * avg)) // если амплитуда достаточно мала или недостаточно велика, то зануляем
        {
            a[i].A = 0;
        }
    }
}

void fill_u_max(vector<dots> &array_fill_from, vector<dots> &array_to_be_filled, unsigned __int64 deltai) 
{
    // заполнение максимумами промежуточного массива u_max и (убрано) подсчет ненулевых элементов в нем, чтобы сформировать окончательный массив
    unsigned __int64 maxindex_from = array_fill_from.size();
    unsigned __int64 maxindex_to_be_filled = array_to_be_filled.size();
    // заполняем массив с максимумами:
    // 1) находим глобальный max
    // 2) записываем его в соотв. массив - array_to_be_filled
    // 3) зануляем элементы вокруг максимума    

    dots min; // глоб минимум
    unsigned __int64 min_index = 0;
    min.A = 0;
    min.t = 0;
    
    for (unsigned __int64 i = 0; i < maxindex_to_be_filled; i++) 
    {
        unsigned __int64 max = findmax(array_fill_from);

        min.A = 0;
        min.t = 0;
        for (unsigned __int64 i = 0; i < maxindex_from; i++) // ищем минимум
        {
            if (array_fill_from[i].A < min.A) 
            {
                min_index = i;
                min = array_fill_from[i];
            }

        }
        
        //cout << "max=" << max << endl;
        if (array_fill_from[max].A > abs(array_fill_from[min_index].A)) // сравниваем по абс. значениям минимум и максимум
        {
            array_to_be_filled[i].t = array_fill_from[max].t;
            array_to_be_filled[i].A = array_fill_from[max].A; // значение за вычетом k*avg       
            del(array_fill_from, max, deltai);
            cout << "fill_u_max: max" << endl;
        }
        else
        {
            array_to_be_filled[i].t = array_fill_from[min_index].t;
            array_to_be_filled[i].A = abs(array_fill_from[min_index].A);            
            del(array_fill_from, min_index, deltai);
            cout << "fill_u_max: min" << endl;
        }       
    }
}

void out(vector<dots> a)  // вывод массива вида a[...].t, a[...].A    
{
    unsigned __int64 maxind = a.size();
    for (unsigned __int64 i = 0; i < maxind; i++) 
    {
        cout << setw(6) << scientific << setprecision(6) << a[i].t << ", " << resetiosflags(ios_base::floatfield);
        cout << setw(6) << setprecision(6) << a[i].A << endl;
    }
    
    cout << resetiosflags(ios_base::floatfield);
}

void file_out(vector<dots> a, filesystem::path p) // запись вектора данных в файл (по сути только для записи ф-и автокорреляции) // еще записывается последний оригинальный сигнал и взвешенный сигнал <- НАДО ПЕРЕДЕЛАТЬ
{   
    // a - вектор данных
    // p - путь файла (вытащили данные из файла p -> записываем обработанные данные в файл с похожим названием)
    unsigned __int64 maxindex = a.size();
    p = p.filename().stem(); // выделяем имя файла без расширения 
    p += "_autocorrelation.txt";
    filesystem::path p1;
    
    p1 = p1 / ".\\computed\\" / p;  // "/" - append

    cout << "file_out:" << p1 << endl;

    ofstream file(p1); // открываем файл для записи

    //cout << fixed << showpoint;
    //cout << setprecision(12);
    file.precision(20);
    //file.setf(ios::fixed);
    file.setf(ios::showpoint);
    for (unsigned __int64 i = 0; i < maxindex; i++)
        file << a[i].t << "," << a[i].A << endl; // записываем
    file.close();
}

void file_out_matrix(vector<dots> max_findmax,             vector<dots> max_autocor,                vector<dots> max_wavelet,               // вектор максимумов, полученный методом поиска макс. и автокор. методом соответствено
                     vector<vector <double> > v_findmax,   vector<vector <double> > v_autocor,      vector<vector <double> > v_wavelet,     // двумерные векторы (матрицы) скоростей
                     vector<vector <double> > at_findmax,  vector<vector <double> > at_autocor,     vector<vector <double> > at_wavelet,    // двумерные векторы (матрицы) коэф. затухания
                     double average_v_findmax,             double average_v_autocor,                double average_v_wavelet,               // средние по таблице
                     double average_attenuation_findmax,   double average_attenuation_autocor,      double average_attenuation_wavelet,
                     double atten_method1_findmax,         double atten_method1_autocor,            double atten_method1_wavelet,           // по методу аппроксимации экспонентой (1)
                     double atten_method2_findmax,         double atten_method2_autocor,            double atten_method2_wavelet,           // по методу аппроксимации экспонентой (2)  
                     filesystem::path path) // путь к файлу, с которым работал алгоритм
{   // вывод в файл

    // p - путь файла (вытащили данные из файла p -> записываем обработанные данные в файл с похожим названием)
    path = path.filename().stem(); // выделяем имя файла без расширения 
    path += "_results.txt";      
    filesystem::path p;
    p = p / ".\\results\\" / path;  // "/" - append
    cout << "file_out_matrix: output file path:" << p << endl;    
    
    ofstream file(p);
    //Проверка успешности открытия файла
    if (file.fail()) {
        cout << "\n Ошибка открытия файла";
        exit(100);
    }

    file.precision(8);    
    file.setf(ios::showpoint);

    file << p.filename().stem() << endl; // выводим имя файла

    file << "***** findmax *****" << endl;
    for (int i = 0; i < max_findmax.size(); i++) // выводим набор максимумов (м. поиска макс.)
        file << "max #" << i << " = " << max_findmax[i].t << "," << max_findmax[i].A << endl;
    file << endl << endl;
    //----
    file << "velocity:" << endl;
    for (int j = 0; j < v_findmax.size(); j++) // проставляем индексы
        file << setw(15) << j; 

    for (int i = 0; i < v_findmax.size(); i++) // выводим скорости (м. поиска макс.)
    {
        file << setw(3) << endl << i; // проставляем индексы
        for (int j = 0; j < v_findmax.size(); j++)
            file << setw(15) << v_findmax[i][j];
    }

    file << endl << endl;
    file << "average velocity: " << endl << average_v_findmax << endl;    
    file << endl << endl;

    //----

    file << "attenuation:" << endl;
    for (int j = 0; j < at_findmax.size(); j++) // проставляем индексы
        file << setw(15) << j;
    
    for (int i = 0; i < at_findmax.size(); i++) // выводим коэф. затух. (м. поиска макс.)
    {
        file << setw(3) << endl << i; // проставляем индексы
        for (int j = 0; j < at_findmax.size(); j++)
            file << setw(15) << at_findmax[i][j];
    }

    file << endl << endl;
    file << "average attenuation coefficient: " << endl << average_attenuation_findmax << endl;
    file << endl << endl;
    file << "attenuation coefficient (method 1): " << endl << atten_method1_findmax << endl;
    file << endl << endl;
    file << "attenuation coefficient (method 2): " << endl << atten_method2_findmax << endl;
    file << endl << endl;

    //----





    file << "***** autocorrelation *****" << endl;
    for (int i = 0; i < max_autocor.size(); i++) // выводим набор максимумов (м. автокор.)
        file << "max #" << i << " = " << max_autocor[i].t << "," << max_autocor[i].A << endl;
    file << endl << endl;
    //----
    file << "velocity:" << endl;
    for (int j = 0; j < v_autocor.size(); j++) // проставляем индексы
        file << setw(15) << j;
    
    for (int i = 0; i < v_autocor.size(); i++) // выводим скорости (м. автокор.)
    {
        file << setw(3) << endl << i; // проставляем индексы
        for (int j = 0; j < v_autocor.size(); j++)
            file << setw(15) << v_autocor[i][j];
    }

    file << endl << endl;
    file << "average velocity: " << endl << average_v_autocor << endl;
    file << endl << endl;


    //----

    file << "attenuation:" << endl;
    for (int j = 0; j < at_autocor.size(); j++) // проставляем индексы
        file << setw(15) << j;

    for (int i = 0; i < at_autocor.size(); i++) // выводим коэф. затух. (м. автокор.)
    {
        file << setw(3) << endl << i; // проставляем индексы
        for (int j = 0; j < at_autocor.size(); j++)
            file << setw(15) << at_autocor[i][j];
    }

    file << endl << endl;
    file << "average attenuation coeffficient: " << endl << average_attenuation_autocor << endl;
    file << endl << endl;
    file << "attenuation coefficient (method 1): " << endl << atten_method1_autocor << endl;
    file << endl << endl;
    file << "attenuation coefficient (method 2): " << endl << atten_method2_autocor << endl;
    file << endl << endl;
    //----



    file << "***** vawelet transform *****" << endl;
    for (int i = 0; i < max_wavelet.size(); i++) // выводим набор максимумов (м. автокор.)
        file << "max #" << i << " = " << max_wavelet[i].t << "," << max_wavelet[i].A << endl;
    file << endl << endl;
    //----
    file << "velocity:" << endl;
    for (int j = 0; j < v_wavelet.size(); j++) // проставляем индексы
        file << setw(15) << j;

    for (int i = 0; i < v_wavelet.size(); i++) // выводим скорости (м. автокор.)
    {
        file << setw(3) << endl << i; // проставляем индексы
        for (int j = 0; j < v_wavelet.size(); j++)
            file << setw(15) << v_wavelet[i][j];
    }

    file << endl << endl;
    file << "average velocity: " << endl << average_v_wavelet << endl;
    file << endl << endl;


    //----

    file << "attenuation:" << endl;
    for (int j = 0; j < at_wavelet.size(); j++) // проставляем индексы
        file << setw(15) << j;

    for (int i = 0; i < at_wavelet.size(); i++) // выводим коэф. затух. (м. автокор.)
    {
        file << setw(3) << endl << i; // проставляем индексы
        for (int j = 0; j < at_wavelet.size(); j++)
            file << setw(15) << at_wavelet[i][j];
    }

    file << endl << endl;
    file << "average attenuation coeffficient: " << endl << average_attenuation_wavelet << endl;
    file << endl << endl;
    file << "attenuation coefficient (method 1): " << endl << atten_method1_wavelet << endl;
    file << endl << endl;
    file << "attenuation coefficient (method 2): " << endl << atten_method2_wavelet << endl;
    file << endl << endl;
    //----



    //----
    //Закрытие файла
    file.close();
}

void check_peak_distance(vector<dots> &array, vector<dots> data, double peak_distance) { // проверка расстояний между пиками
    unsigned __int64 maxindex = array.size();
    // обходим массив и проверяем расстояние между пиками:
    // если расстояние меньше 0,9 * peak_distance , то пик отсеивается;
    //0,9 просто так
     
    // array - массив кандидатов в максимумы
    // data - сигнал

    /*for (int i = 0; i < data.size(); i++)
        cout << data[i].t << ", " << data[i].A << endl;*/


    //unsigned __int64 global_max = findmax(array); // ищем глобальный максимум, чтобы от него отсчитывать расстояние потом

    unsigned __int64 nonzero = 0; // количество ненулевых элементов в array
    for (unsigned __int64 i = 0; i < array.size(); i++) // считаем кол-во ненулевых
        if (array[i].A != 0)
            nonzero++;
    cout << "check_peak_distance: non-zero elements = " << nonzero << endl;

    cout << "check_peak_distance: peak_distance = " << peak_distance << endl;

    unsigned __int64 global_max = floor(nonzero / 2);

    // ищем глобальный максимум в средней части сигнала
    unsigned __int64 l_bound = floor(data.size() * 0.15);
    unsigned __int64 r_bound = floor(data.size() * 0.85);
    cout << "left and right bounds: " << l_bound << ", " << r_bound << endl;

    dots middle_max;
    middle_max.t = 0;
    middle_max.A = 0;

    cout << "check_peak_distance: middle_max (l_bound) = " << middle_max.t << ", " << middle_max.A << endl;
    unsigned __int64 middle_max_index = 0;

    for (int i = l_bound; i < r_bound; i++)
    {
        
        if (data[i].A > middle_max.A)
        {
            middle_max.t = data[i].t;            
            middle_max.A = data[i].A;
            middle_max_index = i;
        }
    }
    cout << "check_peak_distance: middle_max = " << middle_max.t << ", " << middle_max.A << endl;
    //
    // пытаемся найти среди кандидатов в максимумы этот средний максимум, должен найтись по идее
    //
    for (int i = 0; i < array.size(); i++)
    {
        cout << "check_peak_distance: "<< "i = " << i << ", time ratio is "<< abs(array[i].t - middle_max.t) / middle_max.t << endl;
        if ((abs(array[i].t - middle_max.t) / middle_max.t) >= 0 && (abs(array[i].t - middle_max.t) / middle_max.t) < 0.01) // если middle_max примерно равен i-му кандидату в максимум
        {
            global_max = i;
            cout << "check_peak_distance: condition is met; i = " << i << endl;
        }
    }



    /*
    cout << "check_peak_distance: global_max = " << global_max << endl;
    if (((array[global_max + 1].t - array[global_max - 1].t) / 2) < 0.9 * peak_distance || ((array[global_max + 1].t - array[global_max - 1].t) / 2 > 1.1 * peak_distance))
    {
        cout << "check_peak_distance: floor(array.size() / 2) doesn't fit, trying to use floor(array.size() / 3)" << endl;
        global_max = floor(nonzero / 3);
    }
    cout << "check_peak_distance: global_max = " << global_max << endl;
    cout << "u[global_max] = " << array[global_max].t << ", " << array[global_max].A << endl;
    */

    for (unsigned __int64 i = global_max; i < maxindex; i++) { // вперед от максимума
        if (array[i].A == 0) continue;
        if (array[i + 1].t - array[i].t < peak_distance * 0.96)
            array[i + 1].A = 0;
    }

    for (unsigned __int64 i = global_max; i > 0; i--) { // назад от максимума
        if (array[i].A == 0) continue;
        if (array[i].t - array[i - 1].t < peak_distance * 0.96)
            array[i - 1].A = 0;
    }    
}

double z(double d, double a, unsigned __int64 n, unsigned __int64 m) { // параметр z для расчета дифр. поправки; n,m  - номера пиков
    return d / a * sqrt((2 * n - 1) * (2 * m - 1));
}

double c_difr(double c_L, double k_La, double z) { // дифракционная поправка к скорости
    return -c_L * (5.2 / (pow(z, 3 / 2) * k_La * k_La) - 7e-4 * (k_La * k_La - 2200) / (k_La * k_La + 13 * z * z));
}

double calculate_speed(dots maxn, dots maxm, unsigned __int64 n, unsigned __int64 m, double frequency, double specimen_length, double electrode_diameter)
{
    // расчет скорости с учетом дифракционной поправки

    // maxn - первый максимум; n - его номер
    // maxm - второй максимум; m - его номер
    // m>n
    // frequency - частота, нужна для вычисления дифракционной поправки
    // specimen_length - длина (толщина) образца
    // electrode_diameter - диаметр электрода, нужен для выч. дифр. поправки

    if (m < n)
    {
        cout << "calculate_speed: something wrong, m<n" << endl;
        return INFINITY;
    }

    n += 1; // увеличиваем номера, т.к. формулы для поправок написаны так, что нумерация импульсов начинается с 1
    m += 1;

    double speed; // возвращаемое значение скорости
    speed = 2 * specimen_length * (m - n) / (maxm.t - maxn.t); // скорость без дифр. поправки
    speed += c_difr(speed, k_La(frequency, electrode_diameter, speed), z(specimen_length, electrode_diameter, n, m)); // скорость с дифр. поправкой
    //cout << "deltaC_dif=" << c_difr(speed, k_La(frequency, electrode_diameter, speed), z(specimen_length, electrode_diameter, n, m)) << endl;
    //cout << "speed=" << speed;
    
    return speed;
}

complex<double> difrhelp(double k_La, double d, double a, unsigned __int64 n) // вспомогательная ф-я для нахождения дифр. поправки к коэф. затухания
{ 

    // вызывается в ф-и calculate_attenuation, номера импульсов изменяются в ней на 1

    //double xi = k_La / (2 * a) * floor(sqrt((2 * n - 1) * (2 * n - 1) * d * d + 4 * a * a) - (2 * n - 1) * d); // параметр кси для расчета
    //cout << "difrhelp: k=" << k_La / a << endl;    

    double xi = k_La / a / 2 * (sqrt((2 * n - 1) * (2 * n - 1) * d * d + 4 * a * a) - (2 * n - 1) * d); // без floor

    //cout << "difrhelp: xi=" << xi << endl;

    return 1.0 - ((1 - xi * xi / (2 * k_La * k_La)) * _j0(xi) + complex<double>(0, 1) * (1 - xi * xi / (2 * k_La * k_La) + xi / (k_La * k_La)) * _j1(xi)) * exp(-complex<double>(0, 1) * xi);

}

double difr(double k_La, double d, double a, unsigned __int64 n, unsigned __int64 m) // дифр. поправка для коэф. затухания
{ 
    // вызывается в ф-и calculate_attenuation, номера импульсов изменяются в ней на 1

    // k_La - собсно параметр k_La
    // d - толщина образца
    // a - диаметр электрода (?)
    // n, m - номера импульсов
    return 20 * log10(abs(difrhelp(k_La, d, a, n) / difrhelp(k_La, d, a, m)));
}

double calculate_attenuation(double speed, dots maxn, dots maxm, unsigned __int64 n, unsigned __int64 m, double frequency, double specimen_length, double electrode_diameter) // расчет коэф. затухания с учетом дифр. поправки
{ 
    // maxn - первый максимум; n - его номер
    // maxm - второй максимум; m - его номер
    // m>n

    // speed - скорость (для k_La)    
    // frequency - частота, нужна для вычисления дифракционной поправки
    // specimen_length - длина (толщина) образца (d)
    // electrode_diameter - диаметр (?) электрода, нужен для выч. дифр. поправки
    if (m < n)
    {
        cout << "calculate_attenuation: something wrong, m<n" << endl;
        return INFINITY;
    }

    n += 1; // увеличиваем номера, т.к. формулы для поправок написаны так, что нумерация импульсов начинается с 1
    m += 1; 

    double Adifr = difr(k_La(frequency, electrode_diameter, speed), specimen_length, electrode_diameter, n, m); // дифр поправка
    //cout << "calculate_attenuation: difr. correction (dB) = " << Adifr << endl;
    double alpha = (20 * log10(maxn.A / maxm.A) - Adifr) / (2 * specimen_length * (m - n));
    
    //return Adifr;
    return alpha;
}

vector<dots> corresponding_amplitude(vector<dots> a, vector<dots> data, double firstmax) // для нахождения амплитуды, соотв. времени элемента вектора a (чтобы вычислять коэф. затухание автокор. м.)
{
    // a - массив максимумов, для него подбираем амплитуды
    // data - массив данных, в нем ищем соотв. амплитуды
    // firstmax - время "привязки" (надо же от чего-то время отсчитывать)
    vector<dots> result(a.size());

    double t0 = firstmax;
    unsigned __int64 j = 0;
    for (unsigned __int64 i = 0; i < a.size(); i++)
    {
        while (data[j].t < t0 + a[i].t)
            j++;
        result[i].t = a[i].t;
        //result[i].A = (data[j].A + data[j - 1].A + data[j + 1].A + data[j - 2].A + data[j + 2].A) / 5; // среднее по 5 точкам       
        result[i].A = abs(data[j].A);
    }
    return result;
}

vector<dots> autocor_weighing(vector<dots> data, vector<dots> maxima) // взвешиваем набор данных для расчета автокор. ф-и: 
// если элементы массива находятся вблизи одного из максимумов, то вес - 1; если далеко от максимума - меньше 1
// 
{
    double w = 0.0; // вес для далеких от максимумов элементов

    bool near;
    //bool cond;

    double avgd = 0; // средне расстояние между максимумами (по времени)
    for (int i = 0; i < maxima.size() - 1; i++)
        avgd += maxima[i + 1].t - maxima[i].t;
    avgd = avgd / maxima.size();

    cout << endl << endl << "autocor_weighing: avgd = " << avgd << endl << endl;
    double addition = avgd * 0.1875;

    for (unsigned __int64 i = 0; i < data.size(); i++)
    {
        near = false;
        //cond = false;
        for (int j = 0; j < maxima.size(); j++)
        {
            if ((data[i].t > maxima[j].t - addition) && (data[i].t < maxima[j].t + addition)) // max-gap < data < max+gap
            {
                //cond = true;
                near = true;
                //cout << near << endl;
            }
            //else cond = false;
            //near = near || cond;
        }
        //cout << near << endl;
        if (!near)
            data[i].A = data[i].A * w;        

    }
    file_out(data, "weighed_signal"); // записывается с припиской "_autocorrelation", что по смыслу неверно

    return data;
}

double atten_exp_coef(vector<dots> maxima, double velocity, double frequency, double specimen_length, double electrode_diameter) // вычисляем коэф. затухания так: (! и переводим в дБ/м)
// берем максимумы, аппроксимируем их экспонентой по МНК (берем только показатель экспоненты)
// можно применить и другие методы (не МНК)
{
    const double coef = 20 * log10(exp(1)); // коэффициент для перевода Нп -> дБ; умножаем на него
    // обычный МНК для экспоненциальной кривой y=A*exp(Bx)
    double sum_x = 0;
    double sum_lny = 0; 
    double sum_xlny = 0; // sum x*ln(y)
    double sum_x2 = 0;   // sum x^2

    size_t n = maxima.size();

    for (int i = 0; i < n; i++) // вычисляем необходимые суммы
    {
        maxima[i].t = maxima[i].t * velocity; // время*скорость=координата
        sum_x += maxima[i].t;
        sum_lny += log(maxima[i].A);
        sum_xlny += maxima[i].t * log(maxima[i].A);
        sum_x2 += maxima[i].t * maxima[i].t;
    }

    double b = - (n * sum_xlny - sum_x * sum_lny) / (n * sum_x2 - sum_x * sum_x);
    cout << "atten_exp_coef: (without difr correction) alpha (Np) = " << b << endl;
    b *= coef; // переводим из Нп/м в Дб/м
    cout << "atten_exp_coef: (without difr correction) alpha (dB/m) = " << b << endl;
    double Adifr = difr(k_La(frequency, electrode_diameter, velocity), specimen_length, electrode_diameter, 1, n); // дифр поправка; 1 и n - номера импульсов (первый и последний)
    b -= Adifr / (2 * specimen_length * (n - 1)); // вычитаем дифр. поправку, переведенную в дБ/м
    cout << "atten_exp_coef: (with difr correction) alpha (dB/m) = " << b << endl;
    return b;
}

double atten_exp_coef_2(vector<dots> maxima, double velocity, double frequency, double specimen_length, double electrode_diameter) // почти такой же способ, но минимизируется слегка другая функция
// link: https://mathworld.wolfram.com/LeastSquaresFittingExponential.html
{
    const double coef = 20 * log10(exp(1)); // коэффициент для перевода Нп -> дБ; умножаем на него
    size_t n = maxima.size();

    // суммы
    double y = 0;
    double xylny = 0; // x*y*ln(y)
    double xy = 0;
    double ylny = 0;
    double x2y = 0; // x^2 * y

    // можно оптимизировать
    for (int i = 0; i < n; i++) // вычисляем необходимые суммы
    {
        maxima[i].t = maxima[i].t * velocity; // время*скорость=координата

        y += maxima[i].A;
        xylny += maxima[i].t * maxima[i].A * log(maxima[i].A);
        xy += maxima[i].t * maxima[i].A;
        ylny += maxima[i].A * log(maxima[i].A);
        x2y += maxima[i].t * maxima[i].t * maxima[i].A;
    }
    cout << "y=" << y << endl;
    cout << "xylny=" << xylny << endl;
    cout << "xy=" << xy << endl;
    cout << "ylny=" << ylny << endl;
    cout << "x2y=" << x2y << endl;

    double b = - (y * xylny - xy * ylny) / (y * x2y - xy * xy);
    cout << "atten_exp_coef_2: (without difr correction) alpha (Np) = " << b << endl;
    b *= coef; // переводим из Нп/м в Дб/м
    cout << "atten_exp_coef_2: (without difr correction) alpha = " << b << endl;
    double Adifr = difr(k_La(frequency, electrode_diameter, velocity), specimen_length, electrode_diameter, 1, n); // дифр поправка; 1 и n - номера импульсов (первый и последний)
    b -= Adifr / (2 * specimen_length * (n - 1)); // вычитаем дифр. поправку, переведенную в дБ/м
    cout << "atten_exp_coef_2: (with difr correction) alpha = " << b << endl;
    return b;
}

void inversion(vector<dots>& u) // переворачиваем сигнал, если в середине |мин|>|макс| 
{
    // ищем глобальный макс и глобальный мин на участке [0 + 0.15*tlen, tlen - 0.15 * tlen], если |мин|>|макс|, то переворачиваем сигнал отн-но оси времени, т.е. отрицательные будут полож. и наоборот
    // то есть ищем по средней части осциллограммы
    dots globmin = u[floor(u.size() / 2)]; // срединный элемент // просто так
    for (unsigned __int64 i = floor(0.15 * u.size()); i < floor(0.85 * u.size()); i++) // ищем минимум
        if (u[i].A < globmin.A)
            globmin = u[i];
    cout << "inversion: global minimum of [0 + 0.15 * tlen, tlen - 0.15 * tlen] = " << globmin.t << ", " << globmin.A << endl;

    dots globmax = u[floor(u.size() / 2)];
    for (unsigned __int64 i = floor(0.15 * u.size()); i < floor(0.85 * u.size()); i++) // ищем максиум
        if (u[i].A > globmax.A)
            globmax = u[i];
    cout << "inversion: global maximum of [0 + 0.15*tlen, tlen - 0.15 * tlen] = " << globmax.t << ", " << globmax.A << endl;

    if (abs(globmin.A) > abs(globmax.A)) // переворачиваем сигнал
    {
        cout << endl << "inversion: SIGNAL INVERTED" << endl << endl;
        for (unsigned __int64 i = 0; i < u.size(); i++)
        {
            u[i].A *= -1;
        }
    }    
}

vector<dots> maxima(vector<dots> u, string regime, double d, double delta, double f, double a) // получаем максимумы функции
{
    // u - данные для обработки
    // regime - для изменения k дальше   

    vector<dots> u_copy; // копия сигнала
    /*for (unsigned __int64 i = 0; i < u.size(); i++)
        u1.push_back(u[i]);*/
    u_copy = u;
    
    

    unsigned __int64 last = u.size() - 1; // индекс последнего элемента
    double tlen = u[last].t - u[0].t; //длительность осциллограммы по времени

    cout << endl << "tlen=" << tlen << endl;
    //cout << "delta=" << delta << " seconds" << endl;

    

    unsigned __int64 deltai = delta_in_indices(delta, tlen, last); // переводим длительность импульса по времени delta в длину "по индексам"

    //cout << "deltai=" << deltai << " indices" << endl;

    //deltai = deltai_modification(deltai, f); // модифицируем deltai  в зависиомости от частоты: меньше частота - шире пик
    //cout << "modified deltai=" << deltai << " indices" << endl;

    delta = deltai * tlen / (last);
    cout << "(~ in milliseconds) delta=" << delta << endl;

    //------------

    // вычисляем среднюю амплитуду (по абсолютному значению)
    // считаем, не используя ф-ю average, т.к. интересует среднее по абс. значениям
    double avg = 0;
    for (unsigned __int64 i = 0; i < last; i++)
        avg += fabs(u[i].A); //вроде правильно, fabs принимает и выдает double
    avg = avg / (last);

    //cout << "avg=" << avg << endl;

    // вычислили среднее, вывели на экран
    // среднее значение, вычисленное по методу трапеций, не отличается значительно от того, что выше

    // вычислияем "среднее значение" так:
    //      1) вычисляем среднее как выше
    //      2) используя полученное значение, отсеиваем большие по сравнению с этим средним значения
    //      3) получаем новое среднее значение
    //      4) шаг 2
    //      5) и так далее, пока не вступит ограничение на число итераций или на относительную разницу между "новым" и "старым" средними значениями
    // т.о. образуется сходящаяся (?) последовательность значений
    double itavg = iterative_avg(u, 2.8, 0.05, 2);
    // аргументы: u - массив для обработки
    // 5 - коэф. для отсеивания больших значений: отсеиваются значения, в 5 раз больше среднего
    // 0.05 - ограничение на относ. разницу нового и старого средних знач.: если отн. разница меньше 0,05, то итерации прекращаются
    // 10 - ограничение на число итерация: выход из цикла после 10 итерации


    cout << "iterative avg=" << itavg << endl << "avg/itavg=" << avg / itavg << endl;

    avg = itavg;
    //--------------

    double k = 1; //как бы отношение сигнал/шум
    cout << endl << "initial k=" << k << endl;

    k = 0.5*k_modification(k, f); // модифицируем k в зависимости от частоты
    if (regime == "autocorrelation")
        // 2.5 MHz - 1*
        // 10 MHz - 3*
        k = 2*sqrt(k); //  /2 - просто так
        
    cout << "k after modification: k=" << k << endl;


    // k вручную
    // f            k
    // 2.5
    // 5
    // 10           10
    // 16           9
    // 20
    // 32           4
    //cin >> k;

    if (f == 10e6) k = 10;
    if (f == 16e6) k = 9;
    if (f == 32e6) k = 4;


    if (regime == "autocorrelation")
    {
        if (f == 10e6) k = 7;
        if (f == 16e6) k = 11;
    }
    
    cout << "k = " << k <<", k*avg = " << k*avg << endl;
    

    filter_noise(u, k, avg); // отнимаем от всех амлитуд k средних значений, как бы фильтруя сигнал от шума
    cout << "FILTERED" << endl;

    const unsigned __int64 l = 20; // длина вектора с максимумами
    vector<dots> u_max(l); // вектор с максимумами (не окончательный)
    
    //--------------
    // заполняем массив с максимумами:
    // 1) находим глобальный max
    // 2) записываем его в соотв. массив
    // 3) зануляем элементы вокруг максимума
    fill_u_max(u, u_max, deltai); // l - длина массива с максимумами
    out(u_max);
    //--------------

    unsigned __int64 m = 0; // число ненулевых элементов в u_max
    for (unsigned __int64 i = 0; i < l; i++) // подсчитываем число ненулевых элементов в u_max
        if (u_max[i].A > 0)
            m++;
    cout << endl << "m=" << m << endl;

    cout << "u_max=" << endl;
    out(u_max); // значение амплитуды выводится за вычетом k*avg
    cout << endl;

    //------------------------------
    // сортируем пузырьком u_max по времени по (!) убыванию, чтобы узнать расстояние между пиками и затем избавиться от первого пика, если он не годится
    sortd(u_max);
    cout << "descending u_max=" << endl;
    out(u_max);

    // вычисляем среднее расстояние(по времени) между максимумами, чтобы, если случайно 2 макс. "плохие", то это не нарушило бы ничего
    vector<dots> peak_distance(m-1); // <dots> для унификации, по-идее, это лишнее
    for (unsigned __int64 i = 0; i < m - 1; i++) {        
        peak_distance[i].t = u_max[i + 1].t - u_max[i].t;
    }
    double avg_peak_distance = average(peak_distance).t;
    cout << "average peak distance = " << avg_peak_distance << endl;

    //----------    
    sorta(u_max); // сортируем пузырьком u_max по времени по возрастанию, чтобы  потом избавиться от ненужного пика

    cout << "u_max before del=" << endl;
    out(u_max);    

    check_peak_distance(u_max, u_copy, avg_peak_distance);

    cout << "u_max after del=" << endl;
    out(u_max);

    m = 0; //пересчитываем m - количество ненулевых элементов из u_max, чтобы потом сфорировать окончательный массив максимумов
    for (unsigned __int64 i = 0; i < l; i++)
        if (u_max[i].A > 0) m++;
    cout << endl << "m=" << m << endl;
    //-----------------------

    cout << endl << "CHECK #2" << endl << endl;
    for (unsigned __int64 i = 0; i < m - 1; i++) {
        peak_distance[i].t = u_max[i + 1].t - u_max[i].t;
    }
    avg_peak_distance = average(peak_distance).t;
    check_peak_distance(u_max, u_copy, avg_peak_distance);


    cout << "ascending u_max=" << endl;
    out(u_max);

    vector<dots> u_m; // окончательный вектор с данными
    //записываем элементы из u_max в u_m      
    for (unsigned __int64 i = 0; i < l; i++)
        if (u_max[i].A != 0)
            u_m.push_back({ u_max[i].t, u_max[i].A });//+ avg * k }); // прибавляем k*avg, потому что раньше отнимали               
    // записали
    // выводим u_m
    cout << endl << "u_m=" << endl;
    out(u_m);
    cout << endl << endl;

    return u_m;
}

void matr_calc(vector<dots> maxima, vector<vector<double> > &v, vector<vector<double> > &at,  // вычисляем матрицу скоростей и коэф. затухания для набора максимумов и их средние значения
                                    double &avg_v,             double &avg_at           , double f, double d, double a) 
{
    size_t dim = maxima.size();
    for (unsigned __int64 i = 0; i < dim; i++)
        for (unsigned __int64 j = 0; j < dim; j++)
        {
            if (j > i)
            {
                // в формулах нумерация импульсов с 1, встроено в calc_attenuation
                v[i][j] = calculate_speed(maxima[i], maxima[j], i, j, f, d, a); // вычисляем скорость для максимумов i, j
                avg_v += v[i][j];                
                at[i][j] = calculate_attenuation(v[i][j], maxima[i], maxima[j], i, j, f, d, a); // вычисляем коэф. затух. для макс. i, j
                avg_at += at[i][j];
            }
            else
            {
                v[i][j] = 0;
                at[i][j] = 0;
            }
        }

    avg_v = avg_v / (dim * (dim - 1.0) / 2); // ср. арифм. скорости; делитель - количество элементов выше гл. диагонали: N=n*(n-1)/2;
    cout << "matr_calc: avg_v = " << avg_v << endl;

    avg_at = avg_at / (dim * (dim - 1.0) / 2); // ср. арифм. коэф. затухания
    cout << "matr_calc: avg_at = " << avg_at << endl;   
}

complex<double> morlet(double t, double alpha, double omega) { //вейвлет Морле
    //alpha = 1; // пропускная способность // надо попробовать в качестве нее взять коэф затухания    
    //omega = 5 * M_PI;
    //return exp(-t * t / (alpha * alpha)) * exp(complex<double>(0, 1) * omega * t); // проверить скорость вычисления с exp как суммой или произведения
    return exp(-t * t / (alpha * alpha) + complex<double>(0, 1) * omega * t); // вроде быстрее со сложением
}

complex<double> w_morlet_transform(double a, double b, vector<dots> data, double morlet_alpha, double morlet_omega) {
    // a - масштаб (scale)
    // b - сдвиг
    // morlet_alpha и beta - параметры материнского вейвлета Морле
    complex<double> integral = 0;
    for (unsigned __int64 i = 0; i < data.size() - 1; i++) // метод трапеций
    {        
        integral += (conj(morlet(((data[i].t - b) / a), morlet_alpha, morlet_omega)) * data[i].A + conj(morlet(((data[i + 1].t - b) / a), morlet_alpha, morlet_omega)) * data[i + 1].A) * (data[i + 1].t - data[i].t);
    }
    return integral / sqrt(a) / 2.0;
}

double wavelet_max(vector<dots> findmax_maxima, int index, int n, double frequency, double morlet_alpha, vector<dots> data, filesystem::path p)
// поиск максимума при помощи вейвлет-преобразования, чтобы потом использовать его для автокор.
// метода вычисления коэф. затухания и для чего-нибудь еще
// т.е. вместо первого макс. из метода поиска максимумов будет этот максимум

// findmax_maxima - массив максиумов
// index - номер максимума, вокруг которого ищем
// n - число периодов, выдаваемых генератором импульсов
// frequency - несущая частота, Гц 
// n и frequency нужны для определения масштаба a
// data - исходный сигнал вида {t,A} : data <-> u=u(t)
// p - путь к файлу (путь к файлу с исходными данными по идее); нужен чтобы взать имя файла

// !!! масштабный коэф. a определен произвольно, и не совсем понятно как с ним быть

{    
    cout << "wavelet_max: data.front().t = " << data.front().t << ", data.back().t = " << data.back().t << ", data.size() = " << data.size() << endl;

    time_t start, end; // для измерения времени выполнения
    time(&start);    


    // определяем параметры для вейвлет-преобразования W(t,a,b)
    double avg_distance = 0; // средний интервал между максимумами
    for (int i = 0; i < findmax_maxima.size() - 1; i++)
        avg_distance += findmax_maxima[i + 1].t - findmax_maxima[i].t;
    avg_distance = avg_distance / (findmax_maxima.size() - 1); // вычисляем среднее расстояние между максимумами /// <<--- можно ж пробовать брать не среднее расстояние, а туда-сюда на половину длины импульса
    cout << "wavelet_max: avg_distance = " << avg_distance << endl;

    //------------
    double delta = n / frequency;

    //-------------

    // типа выделяем небольшой интервал времени около максимума
    
    //double b_begin = findmax_maxima[index].t - avg_distance / 4;
    double b_begin = findmax_maxima[index].t - 1.1 * delta; // 1,6 вместо 0,5 для запаса
    if (b_begin < 0) // проверяем выход за границы временного интервала сигнала для b_begin
        b_begin = data.front().t; // если <0 получилось, то берем b_begin как начало сигнала
    

    //double b_end   = findmax_maxima[index].t + avg_distance / 4;
    double b_end = findmax_maxima[index].t + 1.1 * delta;
    if (b_end > data.back().t) // проверяем выход за границы временного интервала сигнала для b_end
        b_end = data.back().t; // если получилось больше конечного времени, то берем b_end как конец сигнала
        
    double b_delta = (data[1].t - data.front().t) * 0.5; // шаг -  половина времени дискретизации

    cout << "b before correction:" << endl;
    cout << "b_begin = "    << b_begin << endl << "b_end = " << b_end << endl;
    cout << "b interval = " << b_end - b_begin << " s = "    << (b_end - b_begin) * 1e+9 << " ns" << endl;
    cout << "delta b = "    << b_delta * 1e9   << " ns"      << endl; // в наносекундах

    // делаем так, чтобы начало и конец совпадали с какой-нибудь точкой
    for (int i = 0; i < data.size(); i++)
    {
        if (data[i].t < b_begin) continue;
        b_begin = data[i].t;
        break;
    }

    for (int i = 0; i < data.size(); i++)
    {
        if (data[i].t < b_end) continue;
        b_end = data[i].t;
        break;
    }

    cout << "b after correction:" << endl;
    cout << "b_begin = " << b_begin << endl << "b_end = " << b_end << endl;
    cout << "b interval = " << b_end - b_begin << " s = " << (b_end - b_begin) * 1e+9 << " ns" << endl;
    cout << "delta b = " << b_delta * 1e9 << " ns" << endl; // в наносекундах

    //вычленяем только часть сигнала: ту, куда входит интервал [b_begin, b_end], т.е. выбираем только интересующую часть сигнала
    vector<dots> partial_data;    
    for (unsigned __int64 i = 0; i < data.size(); i++) 
    {
        if (data[i].t < b_begin) continue;
        if (data[i].t > b_end) continue;
        partial_data.push_back(data[i]);
    }
    
    cout << "partial_data.size() = " << partial_data.size() << endl;

    

    // с масштабом непонятно, как быть
    // f=2.5 => a=6e-7..12e-7       a_delta = 4e-9 
    // f=5   => a=3e-7..5e-7        a_delta = 1e-9
    // f=10  => a=1.5e-7..2.75e-7   a_delta = 1e-9 




    double a_begin =0.75*n/2/frequency;
    double a_end = 1.25*n/2/frequency;
    double a_delta = 1e-9;

    cout << "a_begin = " << a_begin << endl << "a_end = " << a_end << endl;
    cout << "a interval = " << a_end - a_begin << endl;
    cout << "delta a = " << a_delta << endl;

    __int64 b_it = ceil((b_end - b_begin) / b_delta); // число итераций по b
    cout << "b: num of itertations = " << b_it << endl;

    __int64 a_it = ceil((a_end - a_begin) / a_delta); // число итераций по a
    cout << "a: num of iterations = " << a_it << endl;

    __int64 ab_it = a_it * b_it;                      // общее число итераций
    cout << "total num of iterations = " << ab_it << endl;


    vector<dwtresult> res (ab_it); // вектор для хранения вейвлет-преобразования
    //cout << "res.size() = " << res.size() << endl << "res.capacity() = "<< res.capacity() << endl;
    

    // настраиваем файл для хранения вейвлет-преобразования (только действительной его части
    // p - путь файла (вытащили данные из файла p -> записываем обработанные данные в файл с похожим названием)    
    p = p.filename().stem(); // выделяем имя файла без расширения 
    cout << "wavelet_max: p=" << p << endl;
    filesystem::path p1;
    p1 += p;
    p1 += "_wavelet_transform_";
    p1 += to_string(index);
    p1 += "_max_re.txt";
    filesystem::path p_re;

    p_re = p_re / ".\\computed" / p1;  // "/" - append

    cout << "wavelet_firstmax: path_re is" << p_re << endl;

    // открываем файл для записи
    ofstream file_re(p_re);  // файл для хранения вейвлет-преобразования (Real part)

    //cout << setprecision(12);
    file_re.precision(10);    
    file_re.setf(ios::showpoint);
    // закончили с файлом
    
    //то же самое делаем для файла для хранения модуля вейвлет-преобразования (просто так, мб пригодится)     
    
    filesystem::path p2;
    p2 += p;
    p2 += "_wavelet_transform_";
    p2 += to_string(index);
    p2 += "_max_abs.txt";
    filesystem::path p_abs;

    p_abs = p_abs / ".\\computed" / p2;  // "/" - append

    cout << "wavelet_firstmax: path_abs is" << p_abs << endl;

    // открываем файл для записи
    ofstream file_abs(p_abs);  // файл для хранения вейвлет-преобразования

    //cout << setprecision(12);
    file_abs.precision(10);
    file_abs.setf(ios::showpoint);
    if (file_abs.fail()) {
        cout << "\n wavelet firstmax: abs File opening error" << endl;
        system("pause");
        exit(101);
    }
    if (file_re.fail()) {
        cout << "\n wavelet firstmax: re File opening error" << endl;
        system("pause");
        exit(102);
    }
    
    // закончили 

    for (__int64 i = 0; i < a_it; i++) // вычисляем вейвлет-преобразование в точке {a,b}; материнский вейвлет - вейвлет Морле (за собственно вычисл. отвечает w_morlet_transform
    {        
        double a = a_begin + i * a_delta;
#pragma omp parallel for schedule(dynamic,100)
        for (__int64 j = 0; j < b_it; j++)
        {            
            double b = b_begin + j * b_delta;            

            //res.push_back({ a_begin + i * a_delta, b_begin + j * b_delta, w_morlet_transform(a,b,partial_data) }); // фигурные скобки для добавления как "точки"
            res[i*b_it + j] = { a_begin + i * a_delta, b_begin + j * b_delta, w_morlet_transform(a, b, partial_data, morlet_alpha, n*M_PI) }; // делаем так, ибо push_back() не является thread-safe методом
            //cout << "a = " << a << ", b = " << b << ", W(a,b)=" << res.back().W << endl;            
            
            if ((i * b_it + j) % 20000 == 0)
            #pragma omp critical (dwt_cout)
            {
                cout << "DWT: iteration # " << i * b_it + j << " of " << ab_it << endl;
                cout << "a = " << a << ", b = " << b << ", W(a,b)=" << res[i * b_it + j].W << endl;
            }
        }
    }


    for (unsigned __int64 i = 0; i < res.size(); i++) // отдельный цикл записи чтобы можно было распараллелить цикл вычислений
    {
        file_re  << res[i].a << " " << res[i].b << " " << res[i].W.real() << endl;
        //file_abs << res[i].a << " " << res[i].b << " " << abs(res[i].W) << endl;
    }

    file_re.close();
    //file_abs.close();
    

    //ищем максимум вейвлет-преобразования (с целью найти максимум в исходном сигнале)
    dwtresult max_re{ res.front().a, res.front().b, res.front().W.real() }; // инициализируем максимум Re(W) как первый элемент res
    dwtresult min_im{ res.front().a, res.front().b, res.front().W.imag() };  // инициализируем минимум Im(W)
    dwtresult max_abs{ res.front().a, res.front().b, abs(res.front().W) };  // инициализируем максимум Abs(W)
    

    for (unsigned __int64 i = 1; i < res.size(); i++) // ищем максиум
    {
        if (res[i].W.real() > max_re.W.real()) // Re
            max_re = res[i];
        if (res[i].W.imag() < min_im.W.imag()) // Im
            min_im = res[i];
        if (abs(res[i].W) > abs(max_abs.W))    // Abs
            max_abs = res[i];
    }
    cout << "wavelet_firstmax: time(b) of real max = "      << max_re.b  << ", a = " << max_re.a  << ", Re(W) = "  << max_re.W.real() << endl;
    cout << "wavelet_firstmax: time(b) of imaginery min = " << min_im.b  << ", a = " << min_im.a  << ", Im(W) = "  << min_im.W.imag() << endl;
    cout << "wavelet_firstmax: time(b) of module max = "    << max_abs.b << ", a = " << max_abs.a << ", abs(W) = " << abs(max_abs.W)  << endl;
    

    time(&end);
    double seconds = difftime(end, start);
    cout << "wavelet_firstmax: runtime = " << seconds / 60 << " minutes" << endl;


    //также пишем в файл некоторые результаты   
    filesystem::path p3;
    p3 += p;
    p3 += "_wavelet_transform_";
    p3 += to_string(index);
    p3 += "_max_re_somedata.txt";
    
    filesystem::path p_data;

    p_data = p_data / ".\\computed" / p3;  // "/" - append

    cout << "wavelet_firstmax: path_data is" << p_data << endl;

    // открываем файл для записи
    ofstream file_data(p_data);  // файл для хранения всяких данных вейвлет-преобразования
    if (file_re.fail()) {
        cout << "\n wavelet firstmax: max_somedata File opening error" << endl;
        system("pause");
        //exit(103);
    }
    // пишем в файл всякие штуки
    file_data << "file name:" << p << endl;
    file_data << "dwt computing time " << seconds << " seconds or " << seconds / 60 << " minutes" << endl;

    file_data << "b_begin = "    << b_begin << endl << "b_end = " << b_end << endl;
    file_data << "b interval = " << b_end - b_begin << " s = "    << (b_end - b_begin) * 1e+9 << " ns" << endl;
    file_data << "delta b = "    << b_delta * 1e9   << " ns"      << endl; // в наносекундах

    file_data << "a_begin = "    << a_begin << endl << "a_end = " << a_end << endl;
    file_data << "a interval = " << a_end - a_begin << endl;
    file_data << "delta a = "    << a_delta << endl;

    file_data << "b: num of itertations = "   << b_it << endl;
    file_data << "a: num of iterations = "    << a_it << endl;
    file_data << "total num of iterations = " << ab_it << endl;

    file_data << "found real maximum: a = "     << max_re.a << ", b = "  << max_re.b  << ", Re(W) = "  << max_re.W.real() << endl;
    file_data << "found imaginary minimum: a = " << min_im.a << ", b = "  << min_im.b  << ", Im(W) = "  << min_im.W.imag() << endl;
    file_data << "found module maximum: a = "   << max_abs.a << ", b = " << max_abs.b << ", abs(W) = " << abs(max_abs.W)  << endl;
    file_data << "difference between findmax_1st_max and wavelet_1st_max (Real) is " << (findmax_maxima[index].t - max_re.b) * 1e9 << " ns" << endl;
    // записали всякую фигню    

    return max_re.b; // возвращаем временной сдвиг - т.е. время максимума
}

int main()
{
    vector<filesystem::path> filenames = get_filenames(".\\waveforms"); // вектор с именами файлов 

    for (int i = 0; i < filenames.size(); i++) // проверяем содержимое вектора с именами файлов
        cout << filenames[i].filename().stem() << filenames[i].extension() << endl; // выводим имя и расширение
    cout << "number of files in folder:" << filenames.size() << endl;
    
    // файл для записи средних значений скорости, коэф. затухания по всем обрабатываемым файлам
    ofstream summary(".\\results\\summary.txt");
    if (summary.fail()) {
        cout << "\n wavelet firstmax: summary file opening error" << endl;
        system("pause");
        //exit(110);
    }

    // заголовки столбцов в итоговом файле c результатами
    summary << setw(40) << " "
        << setw(30) << "avg_C_L_findmax"
        << setw(30) << "avg_atten_findmax"
        << setw(30) << "alpha_findmax"
        << setw(30) << "alpha_findmax_2"

        << setw(30) << "avg_C_L_autocor"
        << setw(30) << "avg_atten_autocor"
        << setw(30) << "alpha_autocor"
        << setw(30) << "alpha_autocor_2"

        << setw(30) << "avg_C_L_wavelet"
        << setw(30) << "avg_atten_wavelet"
        << setw(30) << "alpha_wavelet"
        << setw(30) << "alpha_wavelet_2" << endl;

    // внешние данные
    //-------------------
    double d; // толщина (длина) образца
    cout << "enter the d, meters:" << endl << "d=";
    //cin >> d;
    //d = 0.019948;
    d = 0.01502;
    //d = 0.020454;
    cout << endl;

    
    cout << endl;

    double f; // частота, Гц
    cout << "enter the frequency, Hz" << endl << "f=";
    //cin >> f;
    f = 16e6;
    cout << endl;

    unsigned int n = 8; // число периодов радиоимпульса
    //cin >> n;
    double delta = n / f; // длина радиоимпульса по времени
    

    

    double a; // радиус электрода в метрах
    a = 0.01;
    cout << endl;
    //----------------------

    vector<dots> u; // точки[t,A]

    int filenames_index = 0; // индекс текущего файла
    for (filenames_index = 0; filenames_index < filenames.size(); filenames_index++)    {
        
        
        summary << setw(40) << filenames[filenames_index].filename(); // записываем в файл имена обрабатываемых файлов
        

        cout << endl << "current file: " << filenames[filenames_index] << endl;
        u = fill(filenames[filenames_index], true); // заполнили
        file_out(u, "signal");    

        //inversion(u); // переворачиваем сигнал, если на среднем участке |мин| > |макс|

        cout << filenames[filenames_index] << endl;
        cout << "u[first] = " << u[0].t << ", " << u[0].A << endl;
        cout << "u[last].A = " << u[u.size() - 1].A << endl; // size()-1 <-> back()

        cout << "u_max = " << " # " << findmax(u) << "; A=" << u[findmax(u)].A << endl;

        //----------------------
        // метод поиска максимумов
        // получаем только набор максимумов
        string regime = "findmax";
        vector<dots> findmax_maxima = maxima(u, regime, d, delta, f, a); // получили максимумы из исходных данных
        cout << "***------------***" << endl << "findmax maxima: (" << findmax_maxima.size() << ")" << endl;
        out(findmax_maxima);
        cout << "***------------***" << endl;



        //----------------------
        // автокорреляционный метод
        // получаем только набор максимумов
        regime = "autocorrelation";
        vector<dots> u_weighed = autocor_weighing(u, findmax_maxima); // взвешенный сигнал
        vector<dots> autocor = autocorrelation_calc(u_weighed); //вычисляем автокорреляционную функцию, записываем в вектор autocor
        file_out(autocor, filenames[filenames_index]); // записываем в файл автокорреляционную функцию взвешенного сигнала

        //out(autocor);
        vector<dots> autocor_maxima = maxima(autocor, regime, d, delta, f, a);
        cout << "***------------***" << endl << "autocor maxima: (" << autocor_maxima.size() << ")" << endl;
        out(autocor_maxima);
        cout << "***------------***" << endl;
        
        //------------------



        // вычисляем параметры на основе максимумов (!!!метод поиска максимумов)    
        // делаем это тут, потому что параметры (коэф затухания нужны сразу ниже для вейвлетов)
        cout << endl << "*** computing parameters with findmax method ***" << endl << endl;
        unsigned __int64 dim_findmax = findmax_maxima.size(); // размер вектора findmax_maxima
        //cout << "dim findmax=" << dim_findmax << endl;
        vector<vector<double> > findmax_C_L(dim_findmax, vector<double>(dim_findmax));   // создаем двумерный вектор размером dim_findmax*dim_findmax для расчета скоростей (попарно перебираем все импульсы)
        vector<vector<double> > findmax_atten(dim_findmax, vector<double>(dim_findmax)); // то же самое, только для коэффициента затухания
        //  v_ij   1   2     3     4
        //     1   -   v_21  v_31  v_41
        //     2   -    -    v_32  v_42
        //     3   -    -     -    v_43
        //     4   -    -     -     -

        double avg_C_L_findmax = 0;
        double avg_atten_findmax = 0;

        matr_calc(findmax_maxima, findmax_C_L, findmax_atten, avg_C_L_findmax, avg_atten_findmax, f, d, a);


        double alpha_findmax = atten_exp_coef(findmax_maxima, avg_C_L_findmax, f, d, a);
        cout << "alpha_findmax = " << alpha_findmax << endl;


        double alpha_findmax_2 = atten_exp_coef_2(findmax_maxima, avg_C_L_findmax, f, d, a);
        cout << "alpha_findmax_2 = " << alpha_findmax_2 << endl;

        //-----------------------------------







        // вейвлеты
        //----------------------
        cout << endl << endl;
        vector<dots> wavelet_maxima; //вектор вейвлет-максимумов
        
        //втихушку  определяем коэф затухания (в неперах) для того чтобы вычислить максимумы с помощью вейвлет-преобразования. Коэф. затухания нужен для подстановки его в пропускную способность (по сути коэф. затухания вейвлета) вейвлета Морле
        double wavelet_alpha = alpha_findmax / (20 * log10(exp(1)));
        cout << "alpha for morlet wavelet transform: " << wavelet_alpha << endl;

        for (int i = 0; i < findmax_maxima.size(); i++) // заполняем вектор вейвлет-максимумов 
        {
            wavelet_maxima.push_back({ wavelet_max(findmax_maxima, i, n, f, wavelet_alpha, u, filenames[filenames_index]), 0 }); // ищем время максимумов через вейвлеты, амплитуду полагаем пока равной 0
        }

        if (wavelet_maxima.size() != findmax_maxima.size())
            cout << "Something wrong! Size of wavelet_maxima != size of findmax_maxima" << endl;

        cout << "***------------***" << endl << "wavelet maxima: (" << wavelet_maxima.size() << ")" << endl;
        out(wavelet_maxima);
        cout << "***------------***" << endl;

            

        //------------------

        




        
        
        
        // вычисляем параметры на основе максимумов (!!! через вейвлеты)
        cout << endl << "*** computing parameters with wavelet method ***" << endl << endl;
        unsigned __int64 dim_wavelet = wavelet_maxima.size(); // размер вектора wavelet_maxima

        vector<dots> wavelet_maxima_cor_ampl = corresponding_amplitude(wavelet_maxima, u, 0); // ищем амплитуды, соответствующие найденным точкам по времени при помощи вейвлетов
        cout << "wavelet_maxima_cor_ampl:" << "(" << wavelet_maxima_cor_ampl.size() << ")" << endl;
        cout << "********" << endl;
        out(wavelet_maxima_cor_ampl);
        cout << "********" << endl;

        double avg_C_L_wavelet = 0;
        double avg_atten_wavelet = 0;

        vector<vector<double> > wavelet_C_L(dim_wavelet, vector<double>(dim_wavelet));   // создаем двумерный вектор размером dim_wavelet*dim_wavelet для расчета скоростей (попарно перебираем все импульсы)
        vector<vector<double> > wavelet_atten(dim_wavelet, vector<double>(dim_wavelet)); // то же самое, только для коэффициента затухания



        matr_calc(wavelet_maxima_cor_ampl, wavelet_C_L, wavelet_atten, avg_C_L_wavelet, avg_atten_wavelet, f, d, a);


        double alpha_wavelet = atten_exp_coef(wavelet_maxima_cor_ampl, avg_C_L_wavelet, f, d, a);
        cout << "alpha_wavelet = " << alpha_wavelet << endl;


        double alpha_wavelet_2 = atten_exp_coef_2(wavelet_maxima_cor_ampl, avg_C_L_wavelet, f, d, a);
        cout << "alpha_awavelet_2 = " << alpha_wavelet_2 << endl;
        //----------------------------------------------------
        



        // вычисляем параметры на основе максимумов (!!!автокор. метод) 
        cout << endl << "*** computing parameters with autocorrelation method ***" << endl << endl;
        unsigned __int64 dim_autocor = autocor_maxima.size(); // размер вектора autocor_maxima

        // соотв. амплитуда вычисляется с использованием первого максимума, полученно при помощи метода метода поиска максимумов
        //double wavelet1stmax = wavelet_max(findmax_maxima, 0, u, filenames[i]);
        //cout << "difference between findmax_1st_max and wavelet_1st_max is (dt > 0 -> wavelet1stmax is left):" << (findmax_maxima[0].t - wavelet1stmax) * 1e9 << " ns" << endl;

        vector<dots> autocor_maxima_cor_ampl = corresponding_amplitude(autocor_maxima, u, wavelet_maxima[0].t); // считаем время отн-но 1 максимума из вейвлетов
        //vector<dots> autocor_maxima_cor_ampl = corresponding_amplitude(autocor_maxima, u, findmax_maxima[0].t);        



        cout << "!!!!" << endl;
        out(autocor_maxima_cor_ampl);

        double avg_C_L_autocor = 0;
        double avg_atten_autocor = 0;

        vector<vector<double> > autocor_C_L(dim_autocor, vector<double>(dim_autocor));   // создаем двумерный вектор размером dim_autocor*dim_autocor для расчета скоростей (попарно перебираем все импульсы)
        vector<vector<double> > autocor_atten(dim_autocor, vector<double>(dim_autocor)); // то же самое, только для коэффициента затухания


        matr_calc(autocor_maxima_cor_ampl, autocor_C_L, autocor_atten, avg_C_L_autocor, avg_atten_autocor, f, d, a);


        double alpha_autocor = atten_exp_coef(autocor_maxima_cor_ampl, avg_C_L_autocor, f, d, a);
        cout << "alpha_autocor = " << alpha_autocor << endl;


        double alpha_autocor_2 = atten_exp_coef_2(autocor_maxima_cor_ampl, avg_C_L_autocor, f, d, a);
        cout << "alpha_autocor_2 = " << alpha_autocor_2 << endl;




        file_out_matrix(findmax_maxima,     autocor_maxima_cor_ampl,    wavelet_maxima_cor_ampl,   // векторы максимумов
                        findmax_C_L,        autocor_C_L,                wavelet_C_L,               // 2d-векторы скоростей
                        findmax_atten,      autocor_atten,              wavelet_atten,             // 2d-векторы коэф затухания
                        avg_C_L_findmax,    avg_C_L_autocor,            avg_C_L_wavelet,           // ср. арифм. скоростей 
                        avg_atten_findmax,  avg_atten_autocor,          avg_atten_wavelet,         // ср. арифм. коэф. затухания
                        alpha_findmax,      alpha_autocor,              alpha_wavelet,             // коэф. затухания, полученный как аппроксимация набора максимумов экспонентой 
                        alpha_findmax_2,    alpha_autocor_2,            alpha_wavelet_2,           // то же самое, но аппроксимация слегка другая
                        filenames[filenames_index]);                                               // путь к файлу, с которым работаем

        summary.precision(8);
        summary.setf(ios::showpoint);
        summary
                << setw(30) << avg_C_L_findmax      // ср. скорость (среднее по таблице); м. поиска максимумов
                << setw(30) << avg_atten_findmax    // ср. коэф. затухания (среднее по таблице)
                << setw(30) << alpha_findmax        // коэф. затухания через аппроксимация экспонентой (МНК)
                << setw(30) << alpha_findmax_2      // коэф. затух. также через аппроксимацию экспонентой, но минимизируется слегка другая  ф-я; обычно дает меньшее значение

                << setw(30) << avg_C_L_autocor
                << setw(30) << avg_atten_autocor
                << setw(30) << alpha_autocor
                << setw(30) << alpha_autocor_2

                << setw(30) << avg_C_L_wavelet
                << setw(30) << avg_atten_wavelet
                << setw(30) << alpha_wavelet
                << setw(30) << alpha_wavelet_2 << endl;

    }

    system("pause");
    return 0;
 
}


