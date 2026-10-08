"""
Модуль для работы с биологическими последовательностями в формате FASTA.

Содержит два класса:
    - Seq: хранение и анализ одной последовательности
    - FastaReader: чтение FASTA-файлов (с поддержкой больших файлов через генератор)
"""


import os


class Seq:
    """
    Класс для хранения биологической последовательности (нуклеотидной или белковой).

    Хранит заголовок FASTA-записи и саму последовательность.
    Умеет красиво выводиться в формате FASTA, определять длину
    и тип алфавита (nucleotide / protein / unknown).

    Attributes:
        header (str): Заголовок записи (без символа '>').
        sequence (str): Последовательность в верхнем регистре без пробелов.
    """

    def __init__(self, header, sequence):
        """
        Создаёт объект последовательности.

        Args:
            header (str): Заголовок FASTA-записи.
            sequence (str): Строка с нуклеотидами или аминокислотами.

        Raises:
            TypeError: Если header или sequence не являются строками.
            ValueError: Если header или sequence пустые после очистки.
        """
        if not isinstance(header, str):
            raise TypeError("Header должен быть строкой")
        if not isinstance(sequence, str):
            raise TypeError("Sequence должна быть строкой")
        clean_header = header.strip()
        clean_sequence = "".join(sequence.split()).upper()
        if not clean_header:
            raise ValueError("Заголовок не может быть пустым")
        if not clean_sequence:
            raise ValueError("Введена пустая последовательность")
        self.__header = header.strip()
        self.__sequence = "".join(sequence.split())

    @property
    def header(self):
        """Возвращает заголовок последовательности."""
        return self.__header

    @property
    def sequence(self):
        """Возвращает саму последовательность."""
        return self.__sequence

    def __str__(self):
        """
        Красивый вывод в формате FASTA.

        Заголовок начинается с '>', последовательность разбивается
        на строки по 60 символов.
        """
        lines = [f">{[self.__header]}"]
        for i in range(0, len(self.__sequence), 60):
            lines.append(self.__sequence[i:i+60])
        return "\n".join(lines)


    def __repr__(self):
        """Краткое представление объекта для отладки."""
        return f"Seq(header='{self.__header}', len={len(self.__sequence)})"

    def __len__(self):
        """Возвращает длину последовательности."""
        return len(self.__sequence)

    def alphabet(self):
        """
        Определяет тип алфавита последовательности.

        Returns:
            str: "nucleotide", "protein" или "unknown".
        """
        if not self.__sequence:
            raise ValueError(f"""An empty sequence was introduced!
                             Введена пустая последовательность!""")

        chars = set(self.__sequence)
        nucleotides = set("ACGTUN")
        proteins = set("ACDEFGHIKLMNPQRSTVWY*X")
        if chars.issubset(nucleotides):
            return "nucleotide"
        if chars.issubset(proteins):
            return "protein" 
        return "unknown"

    def is_nucleotide(self):
        """Проверяет, является ли последовательность нуклеотидной."""
        return self.alphabet() == "nucleotide"

    def is_protein(self):
        """Проверяет, является ли последовательность белковой."""
        return self.alphabet() == "protein"



class FastaReader:
    """
    Класс для чтения файлов в формате FASTA.

    Проверяет корректность формата файла и позволяет итерироваться
    по записям, возвращая объекты класса Seq.
    Оптимизирован для работы с большими файлами (использует генератор).

    Attributes:
        path (str): Путь к FASTA-файлу.
    """

    def __init__(self, path):
        """
        Инициализирует читатель FASTA-файла.

        Args:
            path (str): Путь к файлу.

        Raises:
            TypeError: Если path не строка.
            FileNotFoundError: Если файл не существует.
            ValueError: Если путь не является файлом, файл пустой
                        или имеет неверный формат FASTA.
        """

        if not isinstance(path, str):
            raise TypeError("Путь к файлу должен быть строкой")
        if not os.path.exists(path):
            raise FileNotFoundError(f"Файл не найден: {path}")
        if not os.path.isfile(path):
            raise ValueError(f"Указанный путь не является файлом: {path}")
        if os.path.getsize(path) == 0:
            raise ValueError(f"Файл пуст: {path}")

        self.__path = path
        self.__check_format()

    @property
    def path(self):
        """Возвращает путь к файлу."""
        return self.__path

    def __check_format(self):
        """
        Проверяет, что файл начинается с корректной FASTA-записи.

        Raises:
            ValueError: Если файл не начинается с '>' или пустой.
        """

        with open(self.__path, "r", encoding="utf-8") as f:
            for line in f:
                clean_line = line.strip()
                if not clean_line:
                    continue
                if not clean_line.startswith(">"):
                    raise ValueError(f"Неверный формат FASTA: файл должен начинаться с '>', получено: {clean_line[:20]}")
                return
            raise ValueError(f"Файл не содержит данных")

    def read(self):
        """
        Генератор, который читает файл по записям.

        Yields:
            Seq: Объект последовательности для каждой записи в файле.

        Raises:
            ValueError: При обнаружении пустых последовательностей
                        или нарушении формата.
        """
        current_header = None
        current_seq_parts = []

        with open(self.__path, "r", encoding="utf-8") as file:
            for line_num, line in enumerate(file, start=1):
                clean_line = line.strip()

                if not clean_line:
                    continue

                if clean_line.startswith(">"):
                    if current_header is not None:
                        full_seq = "".join(current_seq_parts)
                        if not full_seq:
                            raise ValueError(f"Ошибка в строке {line_num}: пустая последовательность у заголовка {current_header}")
                        yield Seq(current_header, full_seq)

                    current_header = clean_line[1:].strip()
                    current_seq_parts = []
                else:
                    if current_header is None:
                        raise ValueError(f"Ошибка в строке {line_num}: последовательность найдена раньше первого заголовка")
                    current_seq_parts.append(clean_line)

            if current_header is not None:
                full_seq = "".join(current_seq_parts)
                if not full_seq:
                    raise ValueError(f"Пустая последовательность у последнего заголовка: {current_header}")
                yield Seq(current_header, full_seq)

    def __iter__(self):
        """Позволяет использовать объект в цикле for."""
        return self.read()