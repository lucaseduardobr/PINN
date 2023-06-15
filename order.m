% % Ordenar o vetor x em ordem crescente e obter os índices da ordenação
function [x,y] = order(x,y)%,C0
[sorted_x, idx] = sort(x);
 
% % Criar um vetor y com o mesmo tamanho que x e preencher com NaNs
sorted_y = NaN(size(x));
 
% % Atribuir os valores de y que correspondem aos índices ordenados de x
 sorted_y(idx) = y;

end


