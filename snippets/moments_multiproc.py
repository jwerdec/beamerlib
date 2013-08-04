from multiprocessing import Process, Queue, JoinableQueue, freeze_support

def Moments_multiproc(functions, n, num_cpu, limits=(0, np.inf)):
    """
    EXPERIMENTAL
    """
    def worker(in_q, out_q):
        while True:
            indict = in_q.get()
            if indict is None:
                in_q.task_done()
                break
            func = lambda x: functions[indict['funcnum']](x) * x**indict['n']
            indict['moment'] = quad(func, limits[0], limits[1])
            out_q.put(indict)
            in_q.task_done()
        return

    q_in = JoinableQueue()
    q_out = Queue()    
    try:
        iterator = iter(n)
    except TypeError:
        n = [n]
    for funcnum in range(len(functions)):
        for i in n:
            q_in.put({'funcnum': funcnum, 'n': n, 'moment':None})
    for i in range(num_cpu):
        Process(target=worker, args=(q_in, q_out)).start()
    for i in range(num_cpu):
        q_in.put(None)
    q_in.join()
    out = [q_out.get() for _ in range(len(functions)*len(n))]
    result = [{} for func in functions]
    i = 0
    for func in functions:
        for i in n:
            result[i][n] = (item['moment'] for item in out
                            if (item['func']==func) and
                            (n==item['n'])).next()
            i += 1
    return result
