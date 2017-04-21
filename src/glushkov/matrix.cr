struct Matrix(T)
  @array = [] of T
  @row_count : Int32
  @col_count : Int32

  def initialize(@array, @col_count)
    @row_count = @array.size / @col_count
  end

  protected def array
    @array
  end

  def dimensions
    {@row_count, @col_count}
  end

  def same_size?(a : Matrix)
    dim = a.dimensions
    @row_count == dim[0] && @col_count == dim[1]
  end

  def self.convert(m : Matrix)
    new_arr = [] of T
    m.array.each { |e| new_arr << T.new(e) }
    Matrix(T).new(new_arr, m.dimensions[0])
  end

  def self.by_rows(rows)
    col_size = rows[0].size
    raise "wrong size" if rows.any? { |r| r.size != col_size }
    new_arr = [] of T
    rows.each { |r| new_arr.concat(r) }
    new(new_arr, col_size)
  end

  def self.diagonal(e : T, size)
    new_arr = [] of T
    zero = T.zero
    size.times do |i|
      size.times do |j|
        new_arr << (i == j ? e : zero)
      end
    end
    Matrix(T).new(new_arr, size)
  end

  def self.zero(row, col)
    Matrix(T).new(Array.new(row * col, T.zero), col)
  end

  def [](i, j)
    raise "Out of bound" if i >= @row_count || j >= @col_count
    @array[i * @col_count + j]
  end

  def []=(i, j, value : T)
    raise "Out of bound" if i >= @row_count || j >= @col_count
    @array[i * @col_count + j] = value
  end

  def ==(a : Matrix(T))
    @array == a.array
  end

  def +(a : Matrix(T))
    raise "wrong sizes" unless same_size?(a)
    i = -1
    Matrix(T).new(
      a.array.map do |e|
        i += 1
        e + @array[i]
      end,
      a.dimensions[1]
    )
  end

  def -(a : Matrix(T))
    raise "wrong sizes" unless same_size?(a)
    i = -1
    Matrix(T).new(
      a.array.map do |e|
        i += 1
        e - @array[i]
      end,
      a.dimensions[1]
    )
  end

  def *(a : Matrix(T))
    dims = a.dimensions
    raise "wrong sizes" if @row_count != dims[1] || @col_count != dims[0]
    new_arr = [] of T
    @row_count.times do |i1|
      dims[1].times do |j2|
        temp = T.zero
        @col_count.times do |j1|
          temp += self[i1, j1] * a[j1, j2]
        end
        new_arr << temp
      end
    end
    Matrix(T).new(new_arr, dims[1])
  end

  def *(a : T)
    Matrix(T).new(@array.map { |e| e * a }, @col_count)
  end

  def t
    new_arr = [] of T
    @col_count.times do |c|
      @row_count.times do |r|
        new_arr << self[r, c]
      end
    end
    Matrix(T).new(new_arr, @row_count)
  end

  def sub(row, col)
    raise "Out of bounds" if row >= @row_count || col >= @col_count
    new_arr = [] of T
    i = -1
    @row_count.times do |r|
      if r == row
        i += @col_count
        next
      end
      @col_count.times do |c|
        i += 1
        next if c == col
        new_arr << @array[i]
      end
    end
    Matrix(T).new(new_arr, @col_count - 1)
  end

  # minor
  def m(row, col)
    sub(row, col).det
  end

  def cofactor(row, col)
    (-1)**(row + col) * m(row, col)
  end

  def lower_triang
    new_arr = [] of T
    c = -1
    @row_count.times do |i|
      @row_count.times do |j|
        c += 1
        new_arr << (j > i ? T.zero : @array[c])
      end
    end
    Matrix(T).new(new_arr, @row_count)
  end

  def lu_decomposition_crout
    raise "Not square" if @row_count != @col_count
    u = Matrix(Float64).diagonal(1.0, @row_count)
    l = Matrix(Float64).convert(lower_triang)

    # puts self.inspect
    # puts u.inspect
    # puts l.inspect
    # @row_count.times do |j|
    #  next if j == 0
    #  u[0, j] = self[0, j] / l[0, 0]
    # end

    # @row_count.times do |i|
    #  next if i == 0
    #  j = 1
    # #  while j <= i
    sum = 0.0
    @row_count.times do |j|
      i = j
      while i < @row_count
        sum = 0.0
        j.times do |k|
          sum += l[i, k] * u[k, j]
        end
        l[i, j] = self[i, j] - sum
        i += 1
      end

      i = j
      while i < @row_count
        sum = 0.0
        j.times do |k|
          sum += l[j, k] * u[k, i]
        end
        if l[j, j] == 0.0
          raise "det close to zero"
        end
        u[j, i] = (self[j, i] - sum) / l[j, j]
        i += 1
      end
    end

    {u, l}
  end

  def det
    raise "Not square matrix" if @row_count != @col_count
    case @row_count
    when 1
      @array[0]
    when 2
      det_2
    when 3
      det_3
    when 4
      det_4
    else
      u, l = lu_decomposition_crout
      sum1 = u.array[0]
      sum2 = l.array[0]
      @row_count.times do |i|
        next if i == 0
        sum1 *= u[i, i]
        sum2 *= l[i, i]
      end
      sum1 * sum2
    end
  end

  private def det_2
    @array[0] * @array[3] - @array[1] * @array[2]
  end

  private def det_3
    @array[0] * @array[4] * @array[8] +
      @array[1] * @array[5] * @array[6] +
      @array[2] * @array[3] * @array[7] -
      @array[2] * @array[4] * @array[6] -
      @array[1] * @array[3] * @array[8] -
      @array[0] * @array[5] * @array[7]
  end

  private def det_4
    @array[3] * @array[10] * @array[9] * @array[12] - @array[2] * @array[7] * @array[9] * @array[12] -
      @array[3] * @array[5] * @array[10] * @array[12] + @array[1] * @array[7] * @array[10] * @array[12] +
      @array[2] * @array[5] * @array[11] * @array[12] - @array[1] * @array[6] * @array[11] * @array[12] -
      @array[3] * @array[6] * @array[8] * @array[13] + @array[2] * @array[7] * @array[8] * @array[13] +
      @array[3] * @array[4] * @array[10] * @array[13] - @array[0] * @array[7] * @array[10] * @array[13] -
      @array[2] * @array[4] * @array[11] * @array[13] + @array[0] * @array[6] * @array[11] * @array[13] +
      @array[3] * @array[5] * @array[8] * @array[14] - @array[1] * @array[7] * @array[8] * @array[14] -
      @array[3] * @array[4] * @array[9] * @array[14] + @array[0] * @array[7] * @array[9] * @array[14] +
      @array[1] * @array[4] * @array[11] * @array[14] - @array[0] * @array[5] * @array[11] * @array[14] -
      @array[2] * @array[5] * @array[8] * @array[15] + @array[1] * @array[6] * @array[8] * @array[15] +
      @array[2] * @array[4] * @array[9] * @array[15] - @array[0] * @array[6] * @array[9] * @array[15] -
      @array[1] * @array[4] * @array[10] * @array[15] + @array[0] * @array[5] * @array[10] * @array[15]
  end
end
