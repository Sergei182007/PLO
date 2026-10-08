from main import Seq, FastaReader
import sys


def process_file(path):
    print(f"\n{'='*50}")
    print(f"Чтение файла: {path}")
    print('='*50)
   
    try:
        reader = FastaReader(path)
        count = 0
        for record in reader:
            count += 1
            print(f"\n--- Запись #{count} ---")
            print(record)
            print(f"Длина:           {len(record)}")
            print(f"Алфавит:         {record.alphabet()}")
            print(f"is_protein:      {record.is_protein()}")
            print(f"is_nucleotide:   {record.is_nucleotide()}")
       
        if count == 0:
            print("В файле не найдено ни одной записи.")
           
    except FileNotFoundError as e:
        print(f"Файл не найден: {e}")
    except ValueError as e:
        print(f"Ошибка формата/данных: {e}")
    except TypeError as e:
        print(f"Ошибка типа: {e}")
    except Exception as e:
        print(f"Неожиданная ошибка: {e}")


def main():
    if len(sys.argv) > 1:
        for path in sys.argv[1:]:
            process_file(path)
        return

    print("Введите путь к FASTA-файлу")
    user_input = input("> ").strip()

    if user_input:
        paths = user_input.split()
        for path in paths:
            process_file(path)
    else:
        print("Чтение файлов пропущено.")


if __name__ == "__main__":
    main()

