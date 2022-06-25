# Spline curves

compile : g++ main.cpp -I/usr/include/python3.10  -lpython3.10 
OS: LINUX или WSL (Windows Sub-Linux).

Перечень действий по развёртыванию и запуску проекта:
 - Установить пакетным менеджером требуемые модули для Python: 
 		pip3 install matplotlib, numpy, pandas, seaborn;
 - Установить библиотеку python-dev: 
 		sudo apt-get install python-dev (Для Linux Ubuntu);
 - Для тестирования необходимого семейства кривых перейти в каталог на выбор (B-spline, NURBS, X-spline);
 - Перейти в интересующий тест-каталог(см. доступные тесты);
 - Скомпилировать исходный файл утилитой make;
 - Запустить бинарный файл ./a.out с параметрами (см. доступные тесты).



Отрисовка графиков осуществляется с помощью [обёртки matplotlib для С++](https://github.com/lava/matplotlib-cpp) 
Находится в include/matplotlibcpp.h. 
После плоттинга обязательное использование функции PLOT_END().

