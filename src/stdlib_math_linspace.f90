submodule (stdlib_math) stdlib_math_linspace

implicit none

contains


      module procedure linspace_default_1_iint8_iint8

        res = linspace(real(start, kind=dp), real(end, kind=dp), DEFAULT_LINSPACE_LENGTH)

      end procedure linspace_default_1_iint8_iint8
      module procedure linspace_default_1_iint16_iint16

        res = linspace(real(start, kind=dp), real(end, kind=dp), DEFAULT_LINSPACE_LENGTH)

      end procedure linspace_default_1_iint16_iint16
      module procedure linspace_default_1_iint32_iint32

        res = linspace(real(start, kind=dp), real(end, kind=dp), DEFAULT_LINSPACE_LENGTH)

      end procedure linspace_default_1_iint32_iint32
      module procedure linspace_default_1_iint64_iint64

        res = linspace(real(start, kind=dp), real(end, kind=dp), DEFAULT_LINSPACE_LENGTH)

      end procedure linspace_default_1_iint64_iint64

      module procedure linspace_n_1_iint8_iint8

        res = linspace(real(start, kind=dp), real(end, kind=dp), n)

      end procedure linspace_n_1_iint8_iint8
      module procedure linspace_n_1_iint16_iint16

        res = linspace(real(start, kind=dp), real(end, kind=dp), n)

      end procedure linspace_n_1_iint16_iint16
      module procedure linspace_n_1_iint32_iint32

        res = linspace(real(start, kind=dp), real(end, kind=dp), n)

      end procedure linspace_n_1_iint32_iint32
      module procedure linspace_n_1_iint64_iint64

        res = linspace(real(start, kind=dp), real(end, kind=dp), n)

      end procedure linspace_n_1_iint64_iint64

end submodule
