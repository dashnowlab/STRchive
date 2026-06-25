import type { ReactNode } from "react";
import type { JsonData } from "@sagold/json-pointer";
import type { AnnotationData, JsonError, Keyword } from "json-schema-library";
import { cloneElement, Fragment, useMemo } from "react";
import Button from "@/components/Button";
import ComboBox from "@/components/ComboBox";
import { H2 } from "@/components/Heading";
import Help from "@/components/Help";
import NumberBox from "@/components/NumberBox";
import Select from "@/components/Select";
import TextBox from "@/components/TextBox";
import { makeList } from "@/util/format";
import { get, join, remove, set, split } from "@sagold/json-pointer";
import {
  IconArrowDown,
  IconArrowUp,
  IconPlus,
  IconTrash,
} from "@tabler/icons-react";
import clsx from "clsx";
import { compileSchema, draft2020, extendDraft } from "json-schema-library";
import { cloneDeep, isObject, mapValues, range, size, uniq } from "lodash-es";

/** see "format-schema" */
type SchemaNode = {
  $schema?: string;
  $id?: string;
  section?: string;
  title?: string;
  description?: string;
  placeholder?: string;
  examples?: readonly string[];
  "uri-template"?: string;
  type: string | readonly string[];
  pattern?: string;
  minimum?: number;
  maximum?: number;
  less_than?: string;
  greater_than?: string;
  enum?: readonly (string | null)[];
  enum_descriptions?: Record<string, string>;
  in_text_citations?: boolean;
  auto_generated?: boolean;
  hide?: boolean;
  multiline?: boolean;
  combobox?: boolean;
  properties?: Record<string, SchemaNode>;
  uniqueItems?: boolean;
  items?: SchemaNode;
  required?: readonly string[];
  message?: string;
};

type Error = JsonError<
  AnnotationData<{
    message?: string;
    expected?: string;
    pattern?: string;
    maximum?: number;
    minimum?: number;
    minItems?: number;
    maxItems?: number;
  }>
>;

type Props<Schema extends SchemaNode, Data extends JsonData> = {
  schema: Schema;
  sections: string[];
  data: Data;
  onChange: (data: Data) => void;
  children?: ReactNode;
};

/** form automatically generated from json schema */
export default function SchemaForm<
  Schema extends SchemaNode,
  Data extends JsonData,
>({ schema, sections, data, onChange, children }: Props<Schema, Data>) {
  /** compile schema */
  const rootNode = useMemo(
    () => compileSchema(schema, { drafts: [draft] }),
    [schema],
  );

  /** if no initial data, generate blank values */
  data ??= rootNode.getData(null, { extendDefaults: false });

  /** revalidate when data changes */
  const { errors } = useMemo(() => {
    /** @ts-expect-error undocumented way to update data on root node context */
    rootNode.context.data = data; // eslint-disable-line -- see above
    return rootNode.validate(data) as { errors: Error[] };
  }, [rootNode, data]);

  return (
    <>
      {children && <section>{children}</section>}

      {sections.map((section, index) => (
        <section key={index}>
          <H2>{section}</H2>

          <Field<Schema, Data>
            rootNode={rootNode}
            schema={schema}
            section={section}
            node={rootNode}
            data={data}
            setData={onChange}
            errors={errors}
          />
        </section>
      ))}
    </>
  );
}

/** add custom keyword for validating one value less than other value */
const lessThan: Keyword = {
  id: "lessThan",
  keyword: "less_than",
  addValidate: () => true,
  validate: ({ node, data, pointer }) => {
    const thisValue = data;
    /** @ts-expect-error get full data from node context */
    const fullData = node.context.data;
    const otherKey = node.schema.less_than;
    const otherValue = get(fullData, otherKey);

    if (
      typeof thisValue === "number" &&
      typeof otherValue === "number" &&
      thisValue <= otherValue
    )
      return;
    return node.createError("compare-value-error", {
      pointer,
      value: otherValue,
      schema: { ...node.schema, message: `Must be less than ${otherKey}` },
    });
  },
};

/** add custom keyword for validating one value greater than other value */
const greaterThan: Keyword = {
  id: "greaterThan",
  keyword: "greater_than",
  addValidate: () => true,
  validate: ({ node, data, pointer }) => {
    const thisValue = data;
    /** @ts-expect-error undocumented way to get full data from node context */
    const fullData = node.context.data;
    const otherKey = node.schema.greater_than;
    const otherValue = get(fullData, otherKey);
    if (
      typeof thisValue === "number" &&
      typeof otherValue === "number" &&
      thisValue >= otherValue
    )
      return;
    return node.createError("compare-value-error", {
      pointer,
      value: otherValue,
      schema: { ...node.schema, message: `Must be greater than ${otherKey}` },
    });
  },
};

/** https://www.npmjs.com/package/json-schema-library#user-content-extending-a-draft */
const draft = extendDraft(draft2020, {
  /** match all schema */
  $schemaRegEx: ".*",
  /** add custom validation rules */
  keywords: [lessThan, greaterThan],
});

type FieldProps<Schema extends SchemaNode, Data extends JsonData> = {
  rootNode?: ReturnType<typeof compileSchema>;
  schema: Schema;
  section: string;
  node: ReturnType<typeof compileSchema>;
  path?: string;
  data: Data;
  setData: (data: Data) => void;
  errors: Error[];
};

function Field<Schema extends SchemaNode, Data extends JsonData>({
  rootNode,
  schema,
  section,
  node,
  path = "",
  data,
  setData,
  errors,
}: FieldProps<Schema, Data>) {
  /** schema props */
  const {
    title,
    description,
    placeholder,
    examples,
    type,
    enum: _enum,
    enum_descriptions,
    minimum,
    maximum,
    less_than,
    greater_than,
    hide,
    multiline,
    combobox,
  } = node.schema as SchemaNode;

  /** are we at top level of schema */
  const level = split(path).length;

  /** filter out fields that should not be displayed in this section */
  if (level === 1 && (node.schema as SchemaNode).section !== section) return;

  /** explicitly hide field */
  if (hide) return;

  /** normalize type to array */
  const types = [type].flat();

  /** FormData name */
  const name = path;

  /** get nested data value from path */
  const value = get(data, path);

  /** set nested data value from path */
  const onChange = (value: unknown) => {
    let _data = cloneDeep(data);
    _data = set(_data, path, value);
    setData(_data);
  };

  /** help tooltip to show */
  const tooltip = [
    description,
    examples?.length && [
      "Examples:",
      <ul>
        {examples.map((ex) => (
          <li key={ex}>{ex}</li>
        ))}
      </ul>,
    ],
    enum_descriptions &&
      (typeof value === "string" || typeof value === "number") &&
      value in enum_descriptions && [`"${value}": ${enum_descriptions[value]}`],
  ]
    .flat()
    .filter(Boolean)
    .map((line, index) => <div key={index}>{line}</div>);

  /** is the field required to be set */
  const required =
    schema.required?.includes(split(path).pop()!) && !types.includes("null");

  /** full label elements to show */
  const label = (
    <>
      {title ?? path ?? ""}
      {tooltip && <Help>{tooltip}</Help>}
    </>
  );

  /** get validation errors associated with this field */
  const fieldErrors = errors.filter(
    ({ data }) => data.pointer === path.replace(/^#?\/?/, "#/"),
  );
  /** is there an error on this field */
  const isError = !!fieldErrors.length;

  /** error component to show */
  const error = isError ? (
    <span className="col-span-full -mt-2 wrap-break-word text-secondary">
      {fieldErrors
        .map(({ code, data, message }) => {
          const schema = data.schema as SchemaNode;
          const { expected, value, pattern } = data;

          /** prettify error message */
          if (code === "pattern-error")
            return (
              <>
                Must match pattern <code>{pattern}</code>
              </>
            );
          if (code === "enum-error")
            return <>Must be {makeList(schema.enum, "code")}</>;
          if (code === "compare-value-error") return <>{schema.message}</>;
          if (code === "type-error") {
            if (
              [expected].flat().includes("string") &&
              (value === "" || value === null)
            )
              return <>Must not be blank</>;
            else return <>Must be {makeList(expected, "code")}</>;
          }
          if (code === "maximum-error")
            return (
              <>
                Maximum <code>{data.maximum}</code>
              </>
            );
          if (code === "minimum-error")
            return (
              <>
                Minimum <code>{data.minimum}</code>
              </>
            );
          if (code === "unique-items-error") return <>Must not be duplicate</>;
          if (code === "min-items-error")
            return (
              <>
                Must have at least <code>{data.minItems}</code> items
              </>
            );
          if (code === "max-items-error")
            return (
              <>
                Must have at most <code>{data.maxItems}</code> items
              </>
            );
          return message;
        })
        .map((error, index) => (
          <Fragment key={index}>{error}. </Fragment>
        ))}
    </span>
  ) : null;

  /** field control to show */
  let control = (
    <span className="col-span-full -mt-2 wrap-break-word text-secondary">
      missing control
    </span>
  );

  /** ref function */
  const ref = (element: HTMLObjectElement) => {
    if (!element || !("setCustomValidity" in element)) return;
    /** also set native browser form validity to get various built-in features */
    if (isError) element?.setCustomValidity("Error with this field");
    else element?.setCustomValidity("");
  };

  if (types.includes("object"))
    /** object group */
    control = (
      <div
        className={clsx(
          "col-span-full grid w-full grid-cols-[auto_minmax(--spacing(50),1fr)] items-center gap-4 empty:hidden max-md:grid-cols-[auto] [&>label]:contents",
          level > 0 && "rounded-md p-4 shadow-md",
        )}
      >
        {level > 0 && (
          <div className="col-span-full flex items-center gap-2 self-start">
            {label}
          </div>
        )}
        {Object.keys(node.schema.properties).map((key) => {
          /** get child schema node */
          const selection = node.getChildSelection(key);
          if (!Array.isArray(selection) || !selection.length)
            throw Error(`${path} error`);

          return (
            <Field<Schema, Data>
              key={key}
              rootNode={rootNode}
              schema={schema}
              section={section}
              node={{ ...selection[0], schemaAnnotations: [] }}
              path={`${path}/${key}`}
              data={data}
              setData={setData}
              errors={errors}
            />
          );
        })}
      </div>
    );
  else if (types.includes("array")) {
    /** array group */
    const items = size(get(data, path));

    /** get child schema node */
    const selection = node.getChildSelection("items");
    if (!Array.isArray(selection)) throw Error(`${path} error`);

    control = (
      <div className="col-span-full grid w-full grid-flow-row-dense grid-cols-[minmax(--spacing(50),1fr)_auto] items-center gap-4 rounded-md p-4 shadow-md">
        <div className="col-span-full flex items-center gap-2 self-start">
          {label}
        </div>
        {range(items).map((index) => (
          <Fragment key={index}>
            <Field<Schema, Data>
              key={index}
              schema={schema}
              section={section}
              node={{ ...selection[0], schemaAnnotations: [] }}
              path={`${path}/${index}`}
              data={data}
              setData={setData}
              errors={errors}
            />
            <div className="col-start-2 flex flex-wrap items-center *:min-w-10">
              <Button
                design="hollow"
                aria-disabled={index === 0}
                aria-label="Move up"
                onClick={() => {
                  let _data = cloneDeep(data);
                  const aPath = join(path, String(index - 1));
                  const bPath = join(path, String(index));
                  const aValue = get(data, aPath);
                  const bValue = get(data, bPath);
                  _data = set(_data, aPath, bValue);
                  _data = set(_data, bPath, aValue);
                  setData(_data);
                }}
              >
                <IconArrowUp />
              </Button>
              <Button
                design="hollow"
                aria-disabled={index === items - 1}
                aria-label="Move down"
                onClick={() => {
                  let _data = cloneDeep(data);
                  const aPath = join(path, String(index));
                  const bPath = join(path, String(index + 1));
                  const aValue = get(data, aPath);
                  const bValue = get(data, bPath);
                  _data = set(_data, aPath, bValue);
                  _data = set(_data, bPath, aValue);
                  setData(_data);
                }}
              >
                <IconArrowDown />
              </Button>
              <Button
                design="hollow"
                aria-label="Remove"
                onClick={() => {
                  let _data = cloneDeep(data);
                  _data = remove(_data, join(path, String(index)));
                  setData(_data);
                }}
              >
                <IconTrash />
              </Button>
            </div>
          </Fragment>
        ))}

        <Button
          className="col-span-full"
          design="plain"
          onClick={() => {
            /** add new item to array */
            let _data = cloneDeep(data);
            /** get child schema */
            let newItem: unknown = node
              .getNodeChild("items")
              .node?.getData(undefined, {});
            /** modify default value fills */
            if (newItem === "") newItem = null;
            if (isObject(newItem)) newItem = mapValues(newItem, () => null);
            /** insert item */
            _data = set(_data, [...split(path), "[]"], newItem);
            setData(_data);
          }}
        >
          <IconPlus />
          <span>Add</span>
        </Button>
      </div>
    );
  } else if (_enum?.length) {
    /** single select */
    const options = uniq(["", ..._enum.filter(Boolean)]).map((value) => ({
      value: value ?? "",
      label: value,
      tooltip: value !== null ? enum_descriptions?.[value] : "",
    }));
    control = (
      <Select
        options={options}
        value={String(value ?? "")}
        onChange={(value) => onChange(value || null)}
      />
    );
  } else if (types.includes("string"))
    if (combobox && examples) {
      /** text input with autocomplete suggestions */
      control = (
        <ComboBox
          options={examples.map((example) => ({
            value: example,
            label: example,
          }))}
          value={String(value ?? "")}
          onChange={onChange}
        />
      );
    } else
      /** plain text input */
      control = (
        <TextBox
          className="w-full"
          placeholder={
            placeholder ?? (examples?.length ? `Ex: ${examples[0]}` : "")
          }
          multi={multiline}
          value={String(value ?? "")}
          onChange={(value) => onChange(value || null)}
        />
      );
  else if (types.includes("number") || types.includes("integer")) {
    /** min value allowed to be input */
    const min = Math.max(
      ...[minimum, get(data, greater_than ?? "")].filter(
        (value) => typeof value === "number",
      ),
    );
    /** max value allowed to be input */
    const max = Math.min(
      ...[maximum, get(data, less_than ?? "")].filter(
        (value) => typeof value === "number",
      ),
    );
    control = (
      <NumberBox
        step={types.includes("integer") ? 1 : "any"}
        min={Number.isFinite(min) ? min : undefined}
        max={Number.isFinite(max) ? max : undefined}
        value={typeof value === "number" ? value : undefined}
        onChange={onChange}
      />
    );
  }

  return (
    <Fragment>
      {typeof control.type === "function"
        ? cloneElement(
            control,
            /** props that should be on all component controls */
            { ref, name, label, required },
          )
        : control}
      {error}
    </Fragment>
  );
}
