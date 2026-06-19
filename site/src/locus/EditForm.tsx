import type { ReactNode } from "react";
import type { Locus } from "@/data/types";
import { useMemo } from "react";
import { FaXmark } from "react-icons/fa6";
import { LuFeather, LuSend } from "react-icons/lu";
import { createPR } from "@/api/pr";
import Alert from "@/components/Alert";
import Button from "@/components/Button";
import { contactSchema } from "@/components/ContactForm";
import Form from "@/components/Form";
import { H1 } from "@/components/Heading";
import Link from "@/components/Link";
import SchemaForm from "@/components/SchemaForm";
import { repo } from "@/layouts/meta";
import { useQuery } from "@/util/hooks";
import { shortenUrl } from "@/util/string";
import { useLocalStorage } from "@reactuses/core";
import { cloneDeep, isEqual, omitBy, startCase } from "lodash-es";
import loci from "~/STRchive-loci.json";
import _schema from "~/STRchive-loci.schema.json";

/** add extra fields for edit metadata */
const extraProperties = {
  "edit-name": {
    section: "Edit",
    ...contactSchema.name,
    type: ["string", "null"],
    default: "",
  },
  "edit-username": {
    section: "Edit",
    ...contactSchema.username,
    type: ["string", "null"],
    pattern: "^@.+",
    default: "",
  },
  "edit-email": {
    section: "Edit",
    ...contactSchema.email,
    type: ["string", "null"],
    default: "",
  },
  "edit-title": {
    section: "Edit",
    title: "Edit Title",
    description: "Succinct title describing these changes",
    examples: ["Fix mechanism details", "Update disease onset information"],
    type: "string",
    default: null,
  },
  "edit-description": {
    section: "Edit",
    title: "Edit Description",
    description:
      "Summary of changes, justification for changes, uncertainty in literature, or anything else we should know for review. Please be detailed. Provide at least 2-3 sentences.",
    examples: [
      "Currently, the disease mechanism details cite a recently retracted paper doi:123456. This edit corrects the reference and updates...",
    ],
    multiline: true,
    type: "string",
    default: null,
  },
} as const;

/** extra required fields */
const extraRequired = ["edit-title", "edit-description"] as const;

/** add extras */
const schema = {
  ...cloneDeep(_schema),
  properties: { ...extraProperties, ...cloneDeep(_schema.properties) },
  required: [...cloneDeep(_schema.required), ...extraRequired],
};

type Extras = Partial<Record<keyof typeof extraProperties, string | null>>;

type Props = {
  heading?: ReactNode;
  locus?: Locus;
};

/** new/edit locus form */
export default function EditForm({ heading, locus }: Props) {
  /** confirm with user before leaving page */
  // window.onbeforeunload = () => "";

  /** unique storage key for this page and form */
  const storageKey = `edit-locus-${locus?.id ?? "new"}`;

  /** initial form data */
  const initialData = useMemo(
    () =>
      ({
        "edit-title": null,
        "edit-description": null,
        ...cloneDeep(locus),
      }) as Extras & Locus,
    [locus],
  );

  /** form data state */
  const [data, setData] = useLocalStorage(storageKey, initialData);

  /** was data loaded from storage */
  const storageExists = useMemo(() => {
    const fromStorage = window.localStorage.getItem(storageKey);
    /** if saved draft exists and is different from initial data */
    return fromStorage && !isEqual(JSON.parse(fromStorage), initialData);
  }, [storageKey, initialData]);

  /** submission query */
  const {
    query: submit,
    data: response,
    status,
  } = useQuery(async () => {
    /** pr body */
    const name = data?.["edit-name"];
    const username = data?.["edit-username"];
    const email = data?.["edit-email"];
    const description = data?.["edit-description"];
    const body = [{ name, username, email }, { description }]
      .map((group) =>
        Object.entries(group)
          .map(([key, value]) => [
            `**${startCase(key)}**`,
            value?.trim() ? value.trim() : "\\-",
          ])
          .flat()
          .join("\n"),
      )
      .join("\n\n");

    /** branch name from locus id */
    const branch = data?.id;

    /** pr title */
    const title = data?.["edit-title"];

    /** remove edit metadata */
    const newLocus = omitBy(cloneDeep(data), (value, key) =>
      key.startsWith("edit-"),
    ) as Locus;

    /** make complete new clone of loci */
    const newLoci = cloneDeep(loci);

    /** look for existing locus */
    const index = newLoci.findIndex((locus) => locus.id === data?.id);

    /** is new locus vs existing locus */
    const existing = index !== -1;

    /** merge new locus data with existing data */
    if (existing) newLoci[index] = newLocus;
    /** append new locus to end (will be sorted appropriately later) */ else
      newLoci.push(newLocus);

    /** pr files to change */
    const files = [
      {
        path: "data/STRchive-loci.json",
        content: JSON.stringify(newLoci, null, 2),
      },
    ];

    if (!branch) throw Error("branch name missing");
    if (!title) throw Error("title missing");

    return await createPR({
      owner: "dashnowlab",
      repo: "STRchive",
      branch,
      title,
      body,
      files,
      labels: [existing ? "locus-edit" : "locus-new"],
    });
  });

  return (
    <Form onSubmit={submit}>
      <section>
        <H1>
          <LuFeather />
          {heading}
        </H1>

        <Alert type="info">
          Every suggestion is reviewed by our team before inclusion in STRchive.
          Please enter as much accurate information as possible.
        </Alert>

        {storageExists && (
          <div className="flex max-w-full flex-wrap items-center justify-center gap-x-5 gap-y-2">
            Loaded saved draft
            <Button
              design="plain"
              onClick={() => {
                window.localStorage.removeItem(storageKey);
                window.location.reload();
              }}
            >
              Forget
              <FaXmark />
            </Button>
          </div>
        )}
      </section>

      <SchemaForm
        schema={schema}
        data={data}
        onChange={setData}
        sections={[
          "Edit",
          "Overview",
          "Disease",
          "Locus",
          "Alleles",
          "IDs",
          "References",
        ]}
        onSubmit={(data) => submit(data)}
      />

      <section>
        <Alert type={status || "info"}>
          {status === "" && (
            <>
              This will make a <strong>public</strong> pull request on{" "}
              <Link to={repo}>our GitHub</Link>. You'll get a link once it's
              created.
            </>
          )}
          {startCase(status)}{" "}
          {response && (
            <Link to={response.link}>{shortenUrl(response.link)}</Link>
          )}
        </Alert>

        <Button
          type="submit"
          design="bubble"
          aria-disabled={status === "success"}
        >
          <LuSend />
          <span>Submit</span>
        </Button>
      </section>
    </Form>
  );
}
