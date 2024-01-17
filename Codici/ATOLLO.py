import matplotlib.pyplot as plt
import numpy as np
#ciaone

# Dati di esempio per lo scatter plot
x_scatter = np.array([1, 2, 3, 4, 5])
y_scatter = np.array([2, 3, 5, 7, 11])

# Dati di esempio per le funzioni
x_func = np.linspace(0.5, 5.5, 100)  # Creazione di 100 punti tra 0.5 e 5.5
y_func1 = np.sqrt(x_func)
y_func2 = np.log(x_func)

# Creazione del plot scatter
plt.scatter(x_scatter, y_scatter, color='red', label='Scatter Plot')

# Sovrapposizione del grafico delle funzioni
plt.plot(x_func, y_func1, label='Funzione 1')
plt.plot(x_func, y_func2, label='Funzione 2')

# Aggiunta di legenda
plt.legend()

# Mostra il plot
plt.show()