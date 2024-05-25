'''Сохраняет участок файла между start_line и end_line,
если строка после start_line содержит ключевое название гена 
'''

# кортеж ключевых названий генов
brx_names = ('BrxA', 'BrxB', 'BrxC', 'BrxD', 'BrxE', 'BrxF', 'BrxH', 'BrxL', 'BrxP', 'PglW', 'PglX', 'PglZ')

brx_lst = []

start_line = 'HMMER3'
end_line = '//'

# Добавить 'end_line' для гарантированной остановки
with open('padlocdb.hmm', 'a') as file:
    file.write(end_line)

with open('padlocdb.hmm') as file:
    for line in file:
        if line[:6] == start_line:
            line_prev = line
            new_line = file.readline()
            _, name = new_line.split('  ')
            if name[:4] in brx_names:
                brx_lst.append(line_prev)
                brx_lst.append(new_line)
                while not line.startswith(end_line):
                    line = file.readline()
                    brx_lst.append(line)
            else:
                while not line.startswith(end_line):
                    line = file.readline()

with open('brex_select.hmm', mode='w') as file:
    for line in brx_lst:
        file.write(line)
