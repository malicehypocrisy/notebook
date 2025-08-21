# 一、基础

## &#9352;Ndarray 对象

### &#9312;特点

* N维数组对象ndarray，描述相同数据类型的元素集合。
* 以下标0为就开始进行集合元素的索引。

### &#9313;数据类型转换

* astype()函数

## &#9353;数组属性

| 属性             | 说明         |
| ---------------- | ------------ |
| ndarray.ndim     | 秩           |
| ndarray.shape    | 维度         |
| ndarray.size     | 元素的个数   |
| ndarray.dtype    | 元素类型     |
| ndarray.itemsize | 元素的大小   |
| ndarray.flags    | 对象内存信息 |
| ndarray.real     | 元素实部     |
| ndarray.imag     | 元素虚部     |

## &#9354;创建数组

### &#9312;numpy.array()函数

```python
numpy.array(object, dtype = None, copy = True, order = None, subok = False, ndmin = 0)
```

* order：创建数组的形式，C为行方向，F为列方向，A为任意方向（默认）

### &#9313;numpy.zeros()函数

作用：创建默认值为0的数组。

```python
numpy.zeros(shape, dtype = float, order = 'C')
```

* shape:数组形状
* order:有”C“和”F“两个选项，行优先，列优先。

### &#9314;zeros_like()函数

作用：根据传入的ndarray数组的**shape来创建全是零的数组**。

### &#9315;numpy.empty()与empty_like()函数

作用：用于创建空数组，空数据中的值并不为0，而是未初始化的随机值。

```python
numpy.empty(shape, dtype = float, order = 'C')
```

### &#9316;numpy.ones()和ones_like()函数

作用：用于创建所有元素都为1的数组。

### &#9317;numpy.arange()函数

作用：arange函数是python内置函数range函数的数组版本.

```python
numpy.arange(start, stop, step, dtype)
```

### &#9318;numpy.eyes()函数

作用：创建一个N*N的对角矩阵数组，对角线为1，其余为0。

### &#9319;numpy.linspace()函数

作用：创建一个一维数组，该数组是一个等差数列构成。

```python
np.linspace(start, stop, num=50, endpoint=True, retstep=False, dtype=None)
```

* retstep：为TRUE时，生成的数组中会显示间距，反之不显示。

### &#9320;numpy.logspace()函数

```python
np.logspace(start, stop, num=50, endpoint=True, base=10.0, dtype=None)
```

* base：对数log的底数。

## &#9355;从已有数组创建数组

### &#9312;numpy.asarray()函数

```python
numpy.asarray(a, dtype = None, order = None)
```

* a :任意形式的输入参数，可以是列表，列表的元组，元组，元组的元组，元组的列表，多为数组

### &#9313;numpy.frombuffer()函数

作用：实现动态数组

```python
numpy.frombuffer(buffer, dtype = float, count = -1, offset = 0)
```

* count:读取的数据量，默认为-1，读取所有数据。
* offset:读取的其实位置，默认为0。
* buffer是字符串的时候，默认是Unicode类型，转后会在原str前加上b。

### &#9314;numpy.fromiter()函数

作用：从可迭代对象中建立 ndarray 对象，返回一维数组。

```python
numpy.fromiter(iterable, dtype, count=-1)
```

## &#9356;切片与索引

### &#9312;slice()函数

* 设置 start, stop 及 step 参数进行，从原数组中切割出一个新数组

### &#9313;冒号分隔切片

* 不包括停止索引

```python
[start:step:stop]
```

* 切片还可以包括省略号 …，来使选择元组的长度与数组的维度相同。 如果在行位置使用省略号，它将返回包含行中元素的 ndarray。

### &#9314;整数数组索引

* 以坐标的形式索引。
* 借助切片 : 或 … 与索引数组组合
* 布尔索引
* 传入多个索引要使用`np.ix_()函数`

## &#9357;广播

* 如果两个数组a和b形状相同，即满足`a.shape=b.shape`,那么`a*b`就是a与b对应数组位相乘。
* 当维度不同时会触发广播机制（类似R语言的自动补全）。

## &#9358;迭代数组

```python
np.nditer(op, flags=['multi_index'], op_flags=['readwrite'])
```

* `order = 'F'`，即列序优先
* `order = 'C'`，即行序优先

### &#9312;修改数组中元素的值

* 将`op_flags`指定为`read-write`或者`write-only`的模式

### &#9313;广播迭代

* 如果两个数组是可广播的，nditer 组合对象能够同时迭代它们。

## &#9359;数组操作

### &#9312;numpy.reshape()函数

作用：不改变数据的条件下修改形状

# 二、矩阵库

## &#9352;numpy.matlib.empty()函数

作用：返回一个新的矩阵。

```python
numpy.matlib.empty(shape, dtype, order)
```

* `order`:`C`(行序优先)或者`F`(列序优先)

## &#9353;numpy.matlib.zeros()函数

作用：创建一个以 0 填充的矩阵。

## &#9354;numpy.matlib.ones()函数

作用：创建一个以 1 填充的矩阵。

## &#9355;numpy.matlib.eye()函数

作用：返回一个矩阵，对角线元素为1，其他位置为零。

```python
numpy.matlib.eye(n, M,k, dtype)
```

- **n**: 返回矩阵的行数
- **M**: 返回矩阵的列数，默认为 n
- **k**: 对角线的索引
- **dtype**: 数据类型

## &#9356;numpy.matlib.identity()函数

作用：返回给定大小的单位矩阵

## &#9357;numpy.matlib.rand() 函数

作用：创建一个给定大小的矩阵，数据是随机填充的。

# 三、线代