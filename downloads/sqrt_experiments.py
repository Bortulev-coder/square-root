"""Численные эксперименты к РГР «Вычисление квадратного корня».

Python 3.10+. Для графиков: python -m pip install matplotlib
Запуск: python sqrt_experiments.py
Результаты сохраняются в папку sqrt_results рядом с программой.
Все итерации считаются в Decimal с точностью 180 значащих цифр.
sqrt применяется только для эталона, ошибок и проверки результата.
"""
from decimal import Decimal, getcontext
from pathlib import Path
import csv
import json

getcontext().prec = 180
D = Decimal
TOL = D('1e-12')


def heron_step(a, x):
    return (x + a / x) / 2


def bakhshali_step(a, x):
    # Именно исходная формула из задания, не два вызова Герона.
    alpha = (a - x * x) / (2 * x)
    b = x + alpha
    return b - alpha * alpha / (2 * b)


def one_variable(a, x0, steps, method='heron'):
    a, x = D(str(a)), D(str(x0))
    if a <= 0 or x <= 0:
        raise ValueError('Требуются a > 0 и x0 > 0')
    step = {'heron': heron_step, 'bakhshali': bakhshali_step}[method]
    xs = [x]
    for _ in range(steps):
        x = step(a, x)
        xs.append(x)
    return xs


def two_variables(a, steps):
    a = D(str(a))
    if not 0 < a < 3:
        raise ValueError('Требуется 0 < a < 3')
    u, c = a, a - 1  # u обозначает a_n, a — исходное число
    us, cs = [u], [c]
    for _ in range(steps):
        # Обе правые части используют старое значение c.
        u, c = u * (1 - c / 2), c * c * (c - 3) / 4
        us.append(u)
        cs.append(c)
    return us, cs


def relative_errors(xs, a):
    s = D(str(a)).sqrt()
    return [abs(x / s - 1) for x in xs]


def first_hit(xs, a, tol=TOL):
    return next((n for n, e in enumerate(relative_errors(xs, a)) if e <= tol), None)


def residual_stop(x, a, tol=TOL):
    # Достаточный критерий, не использующий sqrt.
    return abs(x * x / a - 1) <= tol


def verify():
    # Допуск существенно выше погрешности 180-разрядной арифметики.
    tiny = D('1e-160')
    for a in map(D, ['0.01', '0.25', '1', '2', '2.9', '2.99']):
        for x0 in map(D, ['0.1', '1', '2', '10']):
            h = one_variable(a, x0, 24)
            b = one_variable(a, x0, 12, 'bakhshali')
            s = a.sqrt()
            assert all(abs(b[n] - h[2*n]) < tiny for n in range(13))
            assert all(x >= s - tiny for x in h[1:])
            assert all(h[n+1] <= h[n] + tiny for n in range(1, 24))
        us, cs = two_variables(a, 35)
        assert all(abs(a*(1+c) - u*u) < tiny for u, c in zip(us, cs))
        assert all(u > 0 for u in us)
        assert all(-1 < c <= 0 for c in cs[1:])
        assert all(cs[n+1] >= cs[n] for n in range(1, 35))
        assert relative_errors(us, a)[-1] < TOL
    return 'Все численные проверки пройдены'


def experiments(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'DejaVu Serif'],
        'mathtext.fontset': 'stix', 'font.size': 11,
        'axes.labelsize': 12, 'axes.titlesize': 12,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.grid': True, 'grid.alpha': .20,
        'figure.dpi': 160, 'savefig.dpi': 240,
    })
    colors = ['#254b70', '#a35038', '#3e7762', '#805b91']
    starts = ['0.1', '1', '2', '10']
    a = D(2)
    s = a.sqrt()
    summary = {'precision': getcontext().prec, 'tolerance': str(TOL)}

    def save(fig, name):
        fig.savefig(out / (name + '.svg'), bbox_inches='tight')
        fig.savefig(out / (name + '.png'), bbox_inches='tight')
        plt.close(fig)

    def errors(ax, xs, radicand, color, label, marker='o'):
        # Не подменяем нулевые ошибки искусственными положительными числами.
        # Точки ниже 1e-100 опускаем: они не нужны для чтения графика.
        pairs = [(n, float(e.log10())) for n,e in enumerate(relative_errors(xs, radicand))
                 if e >= D('1e-100')]
        if pairs:
            nn, yy = zip(*pairs)
            ax.plot(nn, yy, marker=marker, ms=4, color=color, label=label, lw=1.3)
        ax.set_ylabel(r'$\log_{10} E_n$')
        ax.set_xlabel('Номер итерации n')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    for method, prefix, nmax in [('heron','heron',12),('bakhshali','bakhshali',6)]:
        fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.5), layout='constrained')
        for x0, color in zip(starts, colors):
            xs = one_variable(a, x0, nmax, method)
            for ax in axes:
                ax.plot(range(nmax+1), list(map(float,xs)), '-o', color=color,
                        ms=3.5, lw=1.2, label=f'$x_0={x0}$')
                ax.axhline(float(s), color='#222222', ls='--', lw=.7)
                ax.set_xlabel('Номер итерации n')
                ax.set_ylabel(r'$x_n$')
                ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        axes[0].set_title('Общий вид')
        axes[0].legend(fontsize=9)
        axes[1].set_title('Окрестность предела')
        axes[1].set_ylim(1.35,1.65)
        axes[1].set_xlim(1,nmax)
        save(fig, prefix+'_values')
        fig, ax = plt.subplots(figsize=(7.0,3.55),layout='constrained')
        for x0,color in zip(starts,colors):
            errors(ax,one_variable(a,x0,nmax,method),a,color,f'$x_0={x0}$')
        ax.axhline(-12, color='#333333', ls='--', lw=.8, label=r'$E_n=10^{-12}$')
        ax.legend(fontsize=10, ncol=3)
        save(fig,prefix+'_errors')

    values = ['0.01','0.25','2','2.99']
    fig, axes = plt.subplots(2,2,figsize=(7.0,5.5),layout='constrained')
    for ax, aa, color in zip(axes.flat,values,colors):
        us,cs = two_variables(aa,20)
        ax.plot(range(21),list(map(float,us)),'-o',color=color,ms=3,lw=1.2,label=r'$a_n$')
        ax.axhline(float(D(aa).sqrt()),color='#222222',ls='--',lw=.8,label=r'$\sqrt{a}$')
        ax.set_title(f'a = {aa}')
        ax.set_xlabel('Номер итерации n'); ax.set_ylabel(r'$a_n$')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True)); ax.legend(fontsize=10)
    save(fig,'two_values')

    fig,axes=plt.subplots(1,2,figsize=(7.0,3.6),layout='constrained')
    for aa,color in zip(values,colors):
        us,cs=two_variables(aa,23)
        axes[0].plot(range(24),list(map(float,cs)),'-o',color=color,ms=3,lw=1.2,label=f'a={aa}')
        errors(axes[1],us,aa,color,f'a={aa}')
    axes[0].set_xlabel('Номер итерации n'); axes[0].set_ylabel(r'$c_n$')
    axes[0].xaxis.set_major_locator(MaxNLocator(integer=True))
    axes[0].legend(fontsize=9)
    axes[1].axhline(-12,color='#333333',ls='--',lw=.8)
    save(fig,'two_errors')

    counts=[]
    for aa in ['0.01','0.25','1','2','2.9','2.99']:
        h=one_variable(aa,aa,60)
        b=one_variable(aa,aa,30,'bakhshali')
        us,cs=two_variables(aa,60)
        counts.append({'a':aa,'heron':first_hit(h,aa),'bakhshali':first_hit(b,aa),
                       'two_variables':first_hit(us,aa)})
    summary['counts']=counts
    summary['starts']=[{'x0':x0,'heron':first_hit(one_variable(2,x0,30),2),
                        'bakhshali':first_hit(one_variable(2,x0,15,'bakhshali'),2)}
                       for x0 in starts]
    summary['check']=verify()
    (out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    with (out/'iterations.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f); w.writerow(['method','a','x0','n','approximation','relative_error','c_n'])
        for method in ['heron','bakhshali']:
            for x0 in starts:
                xs=one_variable(2,x0,12,method)
                for n,(x,e) in enumerate(zip(xs,relative_errors(xs,2))):
                    w.writerow([method,'2',x0,n,str(x),str(e),''])
        for aa in values:
            us,cs=two_variables(aa,23)
            for n,(x,e,c) in enumerate(zip(us,relative_errors(us,aa),cs)):
                w.writerow(['two_variables',aa,aa,n,str(x),str(e),str(c)])
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__ == '__main__':
    destination=Path(__file__).resolve().parent/'sqrt_results'
    destination.mkdir(exist_ok=True)
    experiments(destination)
